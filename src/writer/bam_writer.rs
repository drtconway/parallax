use std::io::Write;
use std::sync::{Arc, Mutex};

use noodles::bam::io::Writer as BamIoWriter;
use noodles::bgzf::io::Writer as BgzfWriter;
use noodles::sam::{
    Header,
    alignment::{io::Write as _, record_buf::RecordBuf},
};

pub struct BamWriter<W: Write + Send> {
    header: Arc<Header>,
    writer: Mutex<BamIoWriter<BgzfWriter<W>>>,
}

impl<W: Write + Send> BamWriter<W> {
    pub fn new(header: Arc<Header>, writer: W) -> std::io::Result<Self> {
        let writer = noodles::bam::io::Writer::new(writer);
        Ok(BamWriter {
            header,
            writer: Mutex::new(writer),
        })
    }

    pub fn into_inner(self) -> W {
        self.writer.into_inner().unwrap().into_inner().into_inner()
    }
}

impl<W: Write + Send> super::RecordWriter for BamWriter<W> {
    fn write_record(&self, record: &RecordBuf) -> std::io::Result<()> {
        let mut writer = self.writer.lock().unwrap();
        writer.write_alignment_record(&self.header, record)
    }

    fn finish(&self) -> std::io::Result<()> {
        let mut writer = self.writer.lock().unwrap();
        writer.finish(&self.header)
    }
}
