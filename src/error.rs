
pub type Result<T> = std::result::Result<T, ParallaxError>;

#[derive(Debug)]
pub enum ParallaxError {
    IoError(std::io::Error),
    Other(Box<dyn std::error::Error + Send + Sync>),
}

impl std::fmt::Display for ParallaxError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParallaxError::IoError(_) => write!(f, "An IO error occurred."),
            ParallaxError::Other(_) => write!(f, "An error occurred."),
        }
    }
}

impl std::error::Error for ParallaxError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ParallaxError::IoError(e) => Some(e),
            ParallaxError::Other(e) => Some(&**e),
        }
    }
}

impl From<std::io::Error> for ParallaxError {
    fn from(err: std::io::Error) -> Self {
        ParallaxError::IoError(err)
    }
}

impl From<Box<dyn std::error::Error + Send + Sync>> for ParallaxError {
    fn from(err: Box<dyn std::error::Error + Send + Sync>) -> Self {
        ParallaxError::Other(err)
    }
}