#!/usr/bin/env bash
pg_dump ${DBNAME} | xz | gpg --passphrase ${GPG_PASSPHRASE} -c -o dump_${DBNAME}_`date +%d-%m-%Y"_"%H_%M_%S`.xz.gpg

# From docker
#docker exec -t -u postgres pg_dump ${DBNAME} | xz | gpg --passphrase ${GPG_PASSPHRASE} -c -o dump_${DBNAME}_`date +%d-%m-%Y"_"%H_%M_%S`.xz.gpg