# Gaussian Mixture Descriptors Learner (GMDL)

## Install

### Dependencies

Run `sh ./install.sh`. On Ubuntu-based systems run this command as `sudo`.

### App

Run `make` and `./gmdl.app` to run.

## Usage

```
usage: GMDL [options] ... 
options:
  -i, --inline                     pass configurations parameters inline instead of reading json file
  -q, --quiet                      omits logging when classification fails
  -f, --fscore                     outputs f-score instead of confusion matrix
      --no-incremental-learning    disables incremental learning
      --label                      the index of the column that should be considered the label (int [=0])
      --labels                     the labels of the dataset comma-separated (string [=])
      --set                        the key to the set in the config file that should be read (string [=])
  -p, --path                       the path in which to look for datasets (string [=])
      --training                   the name of the traning set (string [=])
      --testing                    the name of the testing set (string [=])
      --config                     the path to the config json file (string [=./config.json])
      --learning_rate              the learning rate used to the SGD on the theta exponents (double [=0.01])
      --momentum                   the momentum used to the SGD on the theta exponents (double [=0.9])
      --tau                        the exponent of the prototype distance used to separate classes (double [=1])
      --beta                       the impact applied to the class that is the most distant from its prototype (double [=-32])
      --omega                      the default description length assumed when there is no clue about the data being assessed (double [=-32])
      --forgetting_factor          the factor by which the samples in the mixture are considered outdated (1 = no forgetting) (double [=1])
      --sigma                      the standard deviation to the noise applied when the covariance matrix is getting close to be singular (double [=1])
  -?, --help                       print this message
```

You can also use a `config.json` with the same keys as in `config.example.json`. Note that you **must** use this json in order to provide the classes being examined.
