from rest_framework import serializers
from dendrogram.models import Profile, Dendrogram


class ProfileSerializer(serializers.ModelSerializer):
    class Meta:
        model = Profile
        fields = ('id', 'batch_id', 'file')


class DendrogramSerializer(serializers.ModelSerializer):
    class Meta:
        model = Dendrogram
        fields = ('id', 'created', 'png_file', 'pdf_file', 'svg_file', 'emf_file', 'newick_file')


class PlotingSerializer(serializers.ModelSerializer):
    class Meta:
        model = Dendrogram
        fields = ('id', )
