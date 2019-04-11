#!/usr/bin/env bash

# Start empty mongodb
docker run --name ${CONTAINER} -d -p 27017:27017 -v ${DATAVOLUME}:/data/db mongo:4.0.5-xenial

# Restore dumped db
docker run --rm --link ${CONTAINER}:mongo -v ${DATAVOLUME}:/backup mongo:4.0.5-xenial \
 bash -c 'mongorestore -d ${DBNAME} --gzip --archive /backup/dump_${DBNAME}_*.archive --host $MONGO_PORT_27017_TCP_ADDR'



# Check data
docker run -it --link some-mongo:mongo --rm mongo:4.0.5-xenial mongo --host mongo vibrio-profiles
