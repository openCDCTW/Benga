import hashlib
import os.path
import datetime
from rest_framework import serializers
from profiling.models import UploadBatch, Sequence, Profile, Dendrogram


class UploadBatchSerializer(serializers.ModelSerializer):
    class Meta:
        model = UploadBatch
        fields = ('id', 'created')


class SequenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = Sequence
        fields = ('id', 'batch_id', 'filename', 'file')

    def create(self, validated_data):
        validated_data["filename"] = hashlib.sha256(validated_data["file"].read()).hexdigest()
        validated_data["file"].name = os.path.join(validated_data["batch_id"],
                                                   validated_data["filename"] + ".fa")
        return Sequence.objects.create(**validated_data)


class ProfileSerializer(serializers.ModelSerializer):
    class Meta:
        model = Profile
        fields = ('id', 'created', 'file', 'occurrence', 'database')

    def create(self, validated_data):
        validated_data["file"].name = "cgmlst-" + str(validated_data["id"].id)[0:8] + ".tsv"
        return Profile.objects.create(**validated_data)


class DendrogramSerializer(serializers.ModelSerializer):
    class Meta:
        model = Dendrogram
        fields = ('id', 'created', 'png_file', 'pdf_file', 'svg_file', 'newick_file')

    def create(self, validated_data):
        now = datetime.datetime.now()
        validated_data["png_file"].name = "cgmlst-" + str(validated_data["id"].id)[0:8] + ".png"
        validated_data["pdf_file"].name = "cgmlst-" + str(validated_data["id"].id)[0:8] + ".pdf"
        validated_data["svg_file"].name = "cgmlst-" + str(validated_data["id"].id)[0:8] + ".svg"
        validated_data["newick_file"].name = "cgmlst-" + str(validated_data["id"].id)[0:8] + ".newick"
        return Dendrogram.objects.create(**validated_data)