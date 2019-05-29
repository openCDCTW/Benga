## Install docker
sudo apt update
sudo apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io
sudo usermod -aG docker $USER


docker pull a504082002/benga:dev



# modify env variables
. envs/web.env
# modify open ip, turn off debug
# check firewall
python3 manage.py makemigrations  # if connection timeout, check the AWS instance security group
python3 manage.py migrate
