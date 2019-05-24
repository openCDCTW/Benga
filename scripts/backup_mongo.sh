#!/usr/bin/env bash
mongodump --db ${DBNAME} --gzip --archive=dump_${DBNAME}_`date +%d-%m-%Y"_"%H_%M_%S`.archive
# mongodump --db ${DBNAME} -o dump_${DBNAME}_`date +%d-%m-%Y"_"%H_%M_%S`

# From docker
#docker run --rm --link some-mongo:mongo -v /data/mongo/backup:/backup mongo:4.0.5-xenial \
# bash -c 'mongodump --db ${DBNAME} --gzip --host $MONGO_PORT_27017_TCP_ADDR --archive /backup/dump_${DBNAME}_`date +%d-%m-%Y"_"%H_%M_%S`.archive'