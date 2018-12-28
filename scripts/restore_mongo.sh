# Start empty mongodb
docker run --name some-mongo -d -p 27017:27017 -v /home/a504082002/Projects/Benga/mongodata:/data/db mongo:4.0.5-xenial

# Restore dumped db
docker run --rm --link some-mongo:mongo -v /home/a504082002/Projects/Benga/vibrio_profiles:/backup mongo:4.0.5-xenial bash -c 'mongorestore /backup/vibrio_profiles -d vibrio-profiles --host $MONGO_PORT_27017_TCP_ADDR'



# Check data
docker run -it --link some-mongo:mongo --rm mongo:4.0.5-xenial mongo --host mongo vibrio-profiles






# Backup DB
docker run \
 --rm \
 --link some-mongo:mongo \
 -v /data/mongo/backup:/backup \
 mongo:4.0.5-xenial \
 bash -c ‘mongodump --out /backup --host $MONGO_PORT_27017_TCP_ADDR’

# Download the dump
scp -r USER@REMOTE:/data/mongo/backup ./backup



# upload file
scp -r ./backup USER@REMOTE:/data/mongo/backup
# Restore DB 
docker run \
 --rm \
 --link some-mongo:mongo \
 -v /data/mongodb/backup:/backup \
 mongo:4.0.5-xenial \
 bash -c ‘mongorestore /backup --host $MONGO_PORT_27017_TCP_ADDR’
