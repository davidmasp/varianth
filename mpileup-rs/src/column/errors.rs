
use thiserror::Error;

#[derive(Error, Debug)]
pub enum NotValidBam {
    #[error("Read name is not set")]
    ReadNameError,
}


