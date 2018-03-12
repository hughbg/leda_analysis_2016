# leda_analysis_2016

Leda data analysis scripts, February 2016


## Install requirements

```
pip install seaborn
pip install lmfit
pip install hickle==2.0.5
pip install bottleneck==0.8.0
```

##Basic usage

Run scripts sequentially, with data file as argument, e.g. `./01_plot_waterfall.py data/outriggers_2017_01_01.h5`.
This is made a bit easier if you set a variable in your terminal like so:

```
> fn=/path/to/data/data.h5
> ./01_plot_waterfall.py $fn
> ./02_plot_spectra.py $fn
```


## Script overview

### 01_plot_waterfall.py

![01-thumb](example_img/01.png =800x)

### 02_plot_spectra.py

![02-thumb](example_img/02.png =800x)

### 02b_plot_spectra_dp.py

![02b-thumb](example_img/02b.png =800x)

### 03_compare_waterfall.py

![03-thumb](example_img/03.png =800x)

### 04_compare_spectra.py

![04-thumb](example_img/04.png =800x)

### 05_plot_residuals.py

![05-thumb](example_img/05.png =800x)

### 06_plot_rfi.py

![06-thumb](example_img/06.png)
![06-thumb](example_img/06_2.png)

### 07_fit_alpha.py

![07-thumb](example_img/07.png)

### 08_plot_fg.py

![08-thumb](example_img/08.png)

### 09_fourier_fitting.py

![09-thumb](example_img/09.png)

### 10_fit_alpha2.py

![10-thumb](example_img/10.png)