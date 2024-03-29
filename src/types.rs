/// Type of the number of points in the metric space.
pub type PointCount = usize;
/// Type of the number of color classes.
pub type ColorCount = usize;
/// Type of the representative intervals [a, b].
pub type Interval = (PointCount,PointCount);
/// Type of the running time return value. Computing time in millisecends.
pub type DurationInSec = f64;

pub(crate) type PointIdx = usize;
pub(crate) type Distance = f32;
pub(crate) type ColorIdx = usize;
pub(crate) type CenterIdx = usize;
