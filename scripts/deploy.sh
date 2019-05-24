#back-end
sudo apt update
sudo apt install git python3-pip virtualenv nginx
git clone https://github.com/openCDCTW/Benga.git
cd Benga/
virtualenv venv -p python3
. venv/bin/activate
pip install --trusted-host pypi.python.org -r requirements.txt --upgrade
#front-end
curl -o- https://raw.githubusercontent.com/creationix/nvm/v0.33.11/install.sh | bash
export NVM_DIR="/root/.nvm"
[ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"
[ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion"
##restart
exit
ssh -i /home/yschen/cdc/cdc-virginia.pem ubuntu@3.87.10.56
cd Benga/
. venv/bin/activate
##
nvm install v10.13.0
npm install
export NODE_OPTIONS=--max_old_space_size=4096
npm run build
. scripts/run_envs.sh  envs/web.env 
. scripts/generate_secret_key.sh
#nginx
pip install uwsgi