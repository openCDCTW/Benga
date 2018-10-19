import json
import os
import shutil
from channels.db import database_sync_to_async
from channels.exceptions import StopConsumer
from channels.generic.websocket import AsyncWebsocketConsumer
from django.conf import settings
from django.core.files import File
from profiling.models import UploadBatch
from profiling.serializers import UploadBatchSerializer, ProfileSerializer, DendrogramSerializer
from src.algorithms import profiling, phylogeny
from src.utils import files


class ProfilingConsumer(AsyncWebsocketConsumer):
    async def connect(self):
        self.batch_id = self.scope['url_route']['kwargs']['batch_id']
        self.channel_id = 'profiling_{}'.format(self.batch_id)

        await self.channel_layer.group_add(
            self.channel_id,
            self.channel_name
        )

        await self.accept()

    async def disconnect(self, close_code):
        await self.channel_layer.group_discard(
            self.channel_id,
            self.channel_name
        )
        raise StopConsumer

    async def receive(self, text_data):
        batch = await self.get_object(self.batch_id)
        serializer = UploadBatchSerializer(batch, data={"id": self.batch_id})
        if serializer.is_valid():
            data_json = json.loads(text_data)
            self.batch = batch

            await self.channel_layer.group_send(
                self.channel_id,
                {
                    "type": "do_profiling",  # processing function
                    "database": data_json["database"],
                    "occr_level": data_json["occr_level"],
                }
            )
        await self.close()

    @database_sync_to_async
    def get_object(self, batch_id):
        try:
            return UploadBatch.objects.get(pk=batch_id)
        except UploadBatch.DoesNotExist:
            return self.close()

    @database_sync_to_async
    def do_profiling(self, event):
        batch = self.batch
        batch_id = str(batch.id)
        database = event["database"]
        occr_level = event["occr_level"]

        input_dir = os.path.join(settings.MEDIA_ROOT, "sequences", batch_id)
        output_dir = os.path.join(settings.MEDIA_ROOT, "temp", batch_id)
        files.create_if_not_exist(output_dir)

        profiling.profiling(output_dir, input_dir, database, occr_level=occr_level, threads=2)

        profile_filename = os.path.join(output_dir, "wgmlst.tsv")

        dendro = phylogeny.Dendrogram()
        dendro.make_tree(profile_filename)
        newick_filename = os.path.join(output_dir, "dendrogram.newick")
        dendro.to_newick(newick_filename)
        pdf_filename = os.path.join(output_dir, "dendrogram.pdf")
        dendro.scipy_tree(pdf_filename)
        svg_filename = os.path.join(output_dir, "dendrogram.svg")
        dendro.scipy_tree(svg_filename)
        png_filename = os.path.join(output_dir, "dendrogram.png")
        dendro.scipy_tree(png_filename)
        # subprocess.call(['libreoffice', '--headless', '--convert-to', 'emf', '--outdir', output_dir, svg_filename])

        profile_data = {"id": batch_id, "file": File(open(profile_filename)),
                        "occurrence": occr_level, "database": database}
        serializer = ProfileSerializer(data=profile_data)
        if serializer.is_valid():
            serializer.save()
        else:
            print(serializer.errors)

        dendrogram_data = {"id": batch_id, "png_file": File(open(png_filename)),
                           "pdf_file": File(open(pdf_filename)), "svg_file": File(open(svg_filename)),
                           "newick_file": File(open(newick_filename))}
        serializer = DendrogramSerializer(data=dendrogram_data)
        if serializer.is_valid():
            serializer.save()
        else:
            print(serializer.errors)

        shutil.rmtree(output_dir)
