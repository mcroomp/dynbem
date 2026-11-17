use pyo3_stub_gen::Result;

fn main() -> Result<()> {
    let stub = _dynbem::stub_info()?;
    stub.generate()?;
    Ok(())
}
