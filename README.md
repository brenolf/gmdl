# Minimum Description Classifier (MDC)

## Install

### Dependencies

Run `sh ./install.sh`. On Ubuntu-based systems run this command as `sudo`.

### App

Run `make` and `./mdc.app` to run.

## Usage

```
usage: MDC [options] ... 
options:
  -i, --inline               pass configurations parameters inline instead of reading json file
  -q, --quiet                omits logging when classification fails
      --label                the index of the column that should be considered the label (int [=0])
      --set                  the key to the set in the config file that should be read (string [=])
      --learning_rate        the learning rate used to the SGD on the theta exponents (double [=0.01])
      --momentum             the momentum used to the SGD on the theta exponents (double [=0.9])
      --delta                the exponent of the prototype distance used to separate classes (double [=1])
      --beta                 the impact applied to the class that is the most distant from its prototype (double [=-32])
      --omega                the default description length assumed when there is no clue about the data being assessed (double [=-32])
      --forgetting_factor    the factor by which the samples in the mixture are considered outdated (1 = no forgetting) (double [=1])
      --sigma                the standard deviation to the noise applied when the covariance matrix is getting close to be singular (double [=1])
  -?, --help                 print this message
```
