# `pwgarch`: An experimental alternative prewhitening scheme for dendro data

## Documentation

See `?pwgarch`.

## Example

```R
library(pwgarch)
data(ca533)
ca_detr <- detrend(ca533, method = "Spline", nyrs = 32)
ca_pw <- pwgarch(ca_detr)
```
