sudo apt update
sudo apt install -y mongodb
# add IP to /etc/mongodb.conf
sudo service mongodb start
sudo systemctl enable mongodb.service
## create user
# use admin
# db.createUser({user: "centrallab",pwd: "5qM5dU5jDf3gVHeP",roles: ["readWrite", "dbAdmin"]})