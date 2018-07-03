#!/bin/sh
echo "Executing R script in the background"
Rscript --vanilla R/forecasting-deforestation-Mada.R > R/forecasting-deforestation-Mada.log 2>&1 &
echo "Check the progress with command 'tail -f R/forecasting-deforestation-Mada.log'"
echo "Check the processor usage with command 'top'"
## End of script