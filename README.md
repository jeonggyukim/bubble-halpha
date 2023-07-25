Below is an example of how you can set up environment. It assumes that you have already installed [miniconda](https://docs.conda.io/en/latest/miniconda.html) or anaconda on your system.

1. Clone the pyathena repo
   ```sh
   git clone https://github.com/jeonggyukim/bubble-halpha.git
   ```
3. Create an environment from the env.yml file
   ```sh
   conda update conda # if you haven't already
   conda env create -f path_to_this_repo/env.yml
   ```
4. Activate the bubble-halpha environment
   ```sh
   conda activate bubble-halpha
   ```
