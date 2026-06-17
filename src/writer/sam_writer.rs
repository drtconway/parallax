use std::io::Write;
use std::sync::{Arc, Mutex};

use noodles::sam::io::Writer as SamIoWriter;
use noodles::sam::{
    Header,
    alignment::{io::Write as _, record_buf::RecordBuf},
};

pub struct SamWriter<W: Write + Send> {
    header: Arc<Header>,
    writer: Mutex<SamIoWriter<W>>,
}

impl<W: Write + Send> SamWriter<W> {
    pub fn new(header: Arc<Header>, writer: W) -> std::io::Result<Self> {
        let writer = noodles::sam::io::Writer::new(writer);
        Ok(SamWriter {
            header,
            writer: Mutex::new(writer),
        })
    }

    pub fn into_inner(self) -> W {
        self.writer.into_inner().unwrap().into_inner()
    }
}

impl<W: Write + Send> super::RecordWriter for SamWriter<W> {
    fn write_record(&self, record: &RecordBuf) -> std::io::Result<()> {
        let mut writer = self.writer.lock().unwrap();
        writer.write_alignment_record(&self.header, record)
    }

    fn finish(&self) -> std::io::Result<()> {
        let mut writer = self.writer.lock().unwrap();
        writer.finish(&self.header)
    }
}
