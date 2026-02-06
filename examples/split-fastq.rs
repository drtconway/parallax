use std::{fs::File, io::Write, path::PathBuf};

use clap::Parser;
use zip::{ZipWriter, write::SimpleFileOptions};

#[derive(Parser, Debug)]
#[command(name = "split-fastq")]
#[command(about = "Split FASTQ files", long_about = None)]
struct Args {
    /// Input FASTQ file (supports gzip/bzip2/xz)
    #[arg(short, long)]
    input: PathBuf,

    /// Output ZIP file
    #[arg(short, long)]
    output: PathBuf,
}

fn main() -> std::io::Result<()> {
    env_logger::builder()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = Args::parse();

    let (decompressed_reader, format) =
        niffler::from_path(&args.input).map_err(|error| std::io::Error::other(error))?;
    if format != niffler::Format::No {
        log::info!("Detected {:?} compression", format);
    }
    let reader = std::io::BufReader::new(decompressed_reader);
    let mut reader = noodles::fastq::io::Reader::new(reader);

    let mut writer = ZipWriter::new(File::create(&args.output)?);
    let options = SimpleFileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated)
        .unix_permissions(0o644);
    for (i, record) in reader.records().enumerate() {
        let record = record?;

        let mut record_writer = noodles::fastq::io::Writer::new(Vec::new());
        record_writer.write_record(&record)?;
        let bytes = record_writer.into_inner();

        let name = String::from_utf8_lossy(record.name());
        let length = record.sequence().len();
        if i % 1_000 == 0 {
            log::info!("Writing read: {}.fastq ({} bases)", name, length);
        }
        let file_name = format!("{}.fastq", name);
        writer.start_file(file_name, options)?;
        writer.write(&bytes)?;
    }
    writer.finish()?;

    Ok(())
}
