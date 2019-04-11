# Restore data to empty db
docker run --name ${CONTAINER} -e POSTGRES_PASSWORD=password -d -p 5432:5432 -v ${DATAVOLUME}:/var/lib/postgresql/data postgres:9.6-alpine

# Restore database
docker exec -i ${CONTAINER} createdb -U postgres ${DBNAME}
gpg --passphrase ${GPG_PASSPHRASE} -o - ${GPGFILE} | xz -d | docker exec -i ${CONTAINER} psql -U postgres ${DBNAME}

docker exec -i ${CONTAINER} createdb -U postgres benga


# Check data
docker run -it --rm --link ${CONTAINER}:postgres postgres:9.6-alpine psql -h postgres -U postgres ${DBNAME}

