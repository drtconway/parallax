use std::sync::Mutex;

use noodles::sam::alignment::RecordBuf;

use crate::writer::RecordWriter;

pub struct SortingWriter<W: RecordWriter> {
    out: W,
    records: Mutex<Vec<RecordBuf>>,
}

impl<W: RecordWriter> SortingWriter<W> {
    pub fn new(out: W) -> Self {
        SortingWriter {
            out,
            records: Mutex::new(Vec::new()),
        }
    }
}

impl<W: RecordWriter> RecordWriter for SortingWriter<W> {
    fn write_record(&self, record: &RecordBuf) -> std::io::Result<()> {
        let mut records = self.records.lock().unwrap();
        records.push(record.clone());
        Ok(())
    }

    fn finish(&self) -> std::io::Result<()> {
        let records = self.records.lock().unwrap();
        let n = records.len();
        let mut order: Vec<usize> = (0..n).collect();

        order.sort_by_key(|&r| {
            let record = &records[r];
            (record.reference_sequence_id(), record.alignment_start())
        });

        for r in order.into_iter() {
            let record = &records[r];
            self.out.write_record(record)?;
        }
        self.out.finish()
    }
}
