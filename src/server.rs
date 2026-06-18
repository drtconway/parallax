use std::io::{BufReader, Cursor};
use std::sync::Arc;

use noodles::sam;
use parallax::{index::Index, reference::InMemoryReference};
use tiny_http::{Header, Method, Request, Response, Server, StatusCode};

use crate::{
    aligner::{Aligner, AlignerBuilder},
    explanatory::ExplanatoryAlignerBuilder,
    writer::{
        RecordWriter,
        bam_writer::BamWriter,
        sam_writer::SamWriter,
        sorting_writer::SortingWriter,
    },
};

const CONTENT_TYPE_BAM: &str = "application/vnd.ga4gh.bam+gzip";
const CONTENT_TYPE_SAM: &str = "text/x-sam";
const CONTENT_TYPE_PLAIN: &str = "text/plain";

enum OutputFormat {
    Bam,
    Sam,
}

pub fn serve(reference: InMemoryReference, index: Arc<dyn Index>, port: u32, secret: Option<String>) {
    let header = Arc::new(build_header(&reference));
    let reference = Arc::new(reference);

    if secret.is_none() {
        log::warn!("no secret set — /quit endpoint is disabled");
    }

    let addr = format!("0.0.0.0:{}", port);
    let server = Server::http(&addr).expect("failed to start HTTP server");
    log::info!("Listening on {}", addr);

    for mut request in server.incoming_requests() {
        let (response, quit) = handle_request(&mut request, &reference, &*index, &header, secret.as_deref());
        if let Err(e) = request.respond(response) {
            log::warn!("failed to send response: {}", e);
        }
        if quit {
            log::info!("shutting down");
            break;
        }
    }
}

fn handle_request(
    request: &mut Request,
    reference: &InMemoryReference,
    index: &dyn Index,
    header: &Arc<sam::Header>,
    secret: Option<&str>,
) -> (Response<Cursor<Vec<u8>>>, bool) {
    let url = request.url().to_string();
    match (request.method(), url.split('?').next().unwrap_or("")) {
        (Method::Post, "/align") => (handle_align(request, reference, index, header), false),
        (Method::Post, "/quit") => handle_quit(&url, secret),
        _ => (text_response(StatusCode(404), "not found"), false),
    }
}

fn handle_quit(url: &str, secret: Option<&str>) -> (Response<Cursor<Vec<u8>>>, bool) {
    let Some(secret) = secret else {
        return (text_response(StatusCode(404), "not found"), false);
    };
    let provided = url
        .split('?')
        .nth(1)
        .and_then(|qs| {
            qs.split('&').find_map(|pair| {
                let (k, v) = pair.split_once('=')?;
                if k == "key" { Some(v) } else { None }
            })
        });
    match provided {
        Some(key) if key == secret => {
            log::info!("quit requested");
            (text_response(StatusCode(200), "ok"), true)
        }
        _ => (text_response(StatusCode(403), "forbidden"), false),
    }
}

fn handle_align(
    request: &mut Request,
    reference: &InMemoryReference,
    index: &dyn Index,
    header: &Arc<sam::Header>,
) -> Response<Cursor<Vec<u8>>> {
    let format = accept_format(request);

    let result = match format {
        OutputFormat::Bam => align_to_bam(reference, index, header, request),
        OutputFormat::Sam => align_to_sam(reference, index, header, request),
    };

    match result {
        Ok(bytes) => {
            let content_type = match format {
                OutputFormat::Bam => CONTENT_TYPE_BAM,
                OutputFormat::Sam => CONTENT_TYPE_SAM,
            };
            bytes_response(bytes, content_type)
        }
        Err(e) => {
            log::error!("alignment failed: {}", e);
            text_response(StatusCode(500), &format!("alignment failed: {}", e))
        }
    }
}

fn accept_format(request: &Request) -> OutputFormat {
    for h in request.headers() {
        if h.field.equiv("accept") {
            if h.value.as_str().contains(CONTENT_TYPE_SAM) {
                return OutputFormat::Sam;
            }
            if h.value.as_str().eq_ignore_ascii_case(CONTENT_TYPE_PLAIN) {
                return OutputFormat::Sam;
            }
        }
    }
    OutputFormat::Bam
}

fn run_aligner<'a>(
    reference: &InMemoryReference,
    index: &dyn Index,
    writer: &dyn RecordWriter,
    request: &mut Request,
) -> std::io::Result<()> {
    let mut aligner = ExplanatoryAlignerBuilder::new(reference, index, writer).build();

    let (decompressed, _format) = niffler::get_reader(Box::new(request.as_reader()))
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    let mut fastq = noodles::fastq::io::Reader::new(BufReader::new(decompressed));

    let mut count = 0;
    for record in fastq.records() {
        let record = record?;
        let name = String::from_utf8_lossy(record.name()).into_owned();
        let seq: &[u8] = record.sequence().as_ref();
        let qual: &[u8] = record.quality_scores().as_ref();
        log::info!("aligning read: {} ({}bp)", name, seq.len());
        aligner.align(&name, seq, qual)?;
        count += 1;
    }
    log::info!("aligned {} read(s), calling aligner.finish()", count);
    aligner.finish()
}

fn align_to_bam<'a>(
    reference: &InMemoryReference,
    index: &dyn Index,
    header: &Arc<sam::Header>,
    request: &mut Request,
) -> std::io::Result<Vec<u8>> {
    let sorting_writer = SortingWriter::new(BamWriter::new(header.clone(), Vec::new())?);
    run_aligner(reference, index, &sorting_writer, request)?;
    Ok(sorting_writer.into_inner().into_inner())
}

fn align_to_sam<'a>(
    reference: &InMemoryReference,
    index: &dyn Index,
    header: &Arc<sam::Header>,
    request: &mut Request,
) -> std::io::Result<Vec<u8>> {
    let sorting_writer = SortingWriter::new(SamWriter::new(header.clone(), Vec::new())?);
    run_aligner(reference, index, &sorting_writer, request)?;
    Ok(sorting_writer.into_inner().into_inner())
}

fn build_header(reference: &InMemoryReference) -> sam::Header {
    use noodles::sam::header::record::value::{Map, map::{self, ReferenceSequence}};
    use noodles::sam::header::record::value::map::header::{sort_order, tag};
    use std::num::NonZero;

    let hd = Map::<map::Header>::builder()
        .insert(tag::SORT_ORDER, sort_order::COORDINATE)
        .build()
        .expect("valid header");

    let mut builder = sam::Header::builder().set_header(hd);
    for (name, len) in reference.chromosomes() {
        if let Some(len) = NonZero::new(len as usize) {
            builder = builder
                .add_reference_sequence(name.to_string(), Map::<ReferenceSequence>::new(len));
        }
    }
    builder.build()
}

fn text_response(status: StatusCode, body: &str) -> Response<Cursor<Vec<u8>>> {
    Response::from_data(body.as_bytes().to_vec())
        .with_status_code(status)
        .with_header(content_type_header("text/plain; charset=utf-8"))
}

fn bytes_response(body: Vec<u8>, content_type: &str) -> Response<Cursor<Vec<u8>>> {
    Response::from_data(body)
        .with_header(content_type_header(content_type))
}

fn content_type_header(value: &str) -> Header {
    Header::from_bytes(b"Content-Type", value.as_bytes()).expect("invalid header")
}
