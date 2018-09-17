import hashlib
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
        validated_data["file"].name = validated_data["filename"] + ".fa"
        return Sequence.objects.create(**validated_data)


class ProfileSerializer(serializers.ModelSerializer):
    class Meta:
        model = Profile
        fields = ('id', 'created', 'file', 'occurrence', 'database')

    def create(self, validated_data):
        validated_data["file"].name = "cgmlst-" + str(validated_data["id"].id) + ".tsv"
        return Profile.objects.create(**validated_data)


class DendrogramSerializer(serializers.ModelSerializer):
    class Meta:
        model = Dendrogram
        fields = ('id', 'created', 'png_file', 'pdf_file', 'svg_file', 'newick_file')

    def create(self, validated_data):
        validated_data["png_file"].name = "cgmlst-" + str(validated_data["id"].id) + ".png"
        validated_data["pdf_file"].name = "cgmlst-" + str(validated_data["id"].id) + ".pdf"
        validated_data["svg_file"].name = "cgmlst-" + str(validated_data["id"].id) + ".svg"
        validated_data["newick_file"].name = "cgmlst-" + str(validated_data["id"].id) + ".newick"
        return Dendrogram.objects.create(**validated_data)