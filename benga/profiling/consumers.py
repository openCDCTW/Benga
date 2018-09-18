from channels.generic.websocket import AsyncWebsocketConsumer
import json


class ProfilingConsumer(AsyncWebsocketConsumer):
    async def connect(self):
        self.room_name = self.scope['url_route']['kwargs']['room_name']
        self.room_group_name = 'chat_%s' % self.room_name

        # Join room group
        await self.channel_layer.group_add(
            self.room_group_name,
            self.channel_name
        )

        await self.accept()

    async def disconnect(self, close_code):
        # Leave room group
        await self.channel_layer.group_discard(
            self.room_group_name,
            self.channel_name
        )

    # Receive message from WebSocket
    async def receive(self, text_data):
        text_data_json = json.loads(text_data)
        message = text_data_json['message']

        # Send message to room group
        await self.channel_layer.group_send(
            self.room_group_name,
            {
                'type': 'chat_message',
                'message': message
            }
        )

    # Receive message from room group
    async def chat_message(self, event):
        message = event['message']

        # Send message to WebSocket
        await self.send(text_data=json.dumps({
            'message': message
        }))

    # async def do_profiling(self, batch_id, database, occr_level):
    #     input_dir = os.path.join(INDIR, batch_id)
    #     files.create_if_not_exist(OUTDIR)
    #     output_dir = os.path.join(OUTDIR, batch_id)
    #     files.create_if_not_exist(output_dir)
    #     profiling.profiling(output_dir, input_dir, database, occr_level=occr_level, threads=2)
    #     profile_created = datetime.datetime.now()
    #
    #     with open(os.path.join(output_dir, "namemap.json"), "r") as file:
    #         names = json.loads(file.read())
    #     profile_filename = os.path.join(output_dir,
    #                                     "cgMLST_{}_{}_{}.tsv".format(database, occr_level, batch_id[0:8]))
    #     os.rename(os.path.join(output_dir, "wgmlst.tsv"), profile_filename)
    #     dendro = phylogeny.Dendrogram()
    #     dendro.make_tree(profile_filename, names)
    #     dendro_created = datetime.datetime.now()
    #     newick_filename = os.path.join(output_dir, "dendrogram_{}.newick".format(batch_id[0:8]))
    #     dendro.to_newick(newick_filename)
    #     pdf_filename = os.path.join(output_dir, "dendrogram_{}.pdf".format(batch_id[0:8]))
    #     dendro.scipy_tree(pdf_filename)
    #     svg_filename = os.path.join(output_dir, "dendrogram_{}.svg".format(batch_id[0:8]))
    #     dendro.scipy_tree(svg_filename)
    #     png_filename = os.path.join(output_dir, "dendrogram_{}.png".format(batch_id[0:8]))
    #     dendro.scipy_tree(png_filename)
    #
    #     sql = "INSERT INTO profile (id,created,file,occurrence,database)" \
    #           "VALUES(:id,:created,:file,:occr,:db);"
    #     data = {"id": batch_id, "created": profile_created, "file": profile_filename,
    #             "occr": occr_level, "db": database}
    #     db.to_sql(sql, data)
    #
    #     sql = "INSERT INTO dendrogram (id,created,png_file,pdf_file,svg_file,newick_file)" \
    #           "VALUES(:id,:created,:png,:pdf,:svg,:new);"
    #     data = {"id": batch_id, "created": dendro_created, "png": png_filename, "pdf": pdf_filename,
    #             "svg": svg_filename, "new": newick_filename}
    #     db.to_sql(sql, data)