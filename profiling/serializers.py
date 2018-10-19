from rest_framework import serializers
from profiling.models import UploadBatch, Sequence, Profile, Dendrogram


class UploadBatchSerializer(serializers.ModelSerializer):
    class Meta:
        model = UploadBatch
        fields = ('id', 'created')


class SequenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = Sequence
        fields = ('id', 'batch_id', 'file')


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
