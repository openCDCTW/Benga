# Restore data to empty db
docker run --name some-postgres -e POSTGRES_PASSWORD=password -d -p 5432:5432 -v /home/a504082002/Projects/Benga/pgdata:/var/lib/postgresql/data postgres:9.6-alpine

# Restore database
docker exec -i some-postgres createdb -U postgres Vibrio_cholerae
gunzip -c 20181227_vibrio.gz | docker exec -i some-postgres psql -U postgres Vibrio_cholerae

docker exec -i some-postgres createdb -U postgres benga


# Check data
docker run -it --rm --link some-postgres:postgres postgres:9.6-alpine psql -h postgres -U postgres Vibrio_cholerae






# Backup
docker exec -t -u postgres your-db-container pg_dumpall -c > dump_`date +%d-%m-%Y"_"%H_%M_%S`.sql

# Restore
cat your_dump.sql | docker exec -i your-db-container psql -U postgres
