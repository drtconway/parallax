use std::io::Write;

pub trait DumpItem: Sized {
    type HeaderInfo;

    fn write_header(header_info: &Self::HeaderInfo, writer: &mut impl Write);

    fn write(&self, writer: &mut impl Write);
}
