from django.http import Http404
from rest_framework import generics, status
from rest_framework.response import Response
from rest_framework.views import APIView

from dendrogram.models import Batch, Profile, Dendrogram
from dendrogram.serializers import BatchSerializer, ProfileSerializer,\
    DendrogramSerializer, PlotingSerializer
from dendrogram.tasks import plot_dendrogram


class BatchList(generics.ListCreateAPIView):
    queryset = Batch.objects.all()
    serializer_class = BatchSerializer


class BatchDetail(APIView):
    def get_object(self, pk):
        try:
            return Batch.objects.get(pk=pk)
        except Batch.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        batch = self.get_object(pk)
        serializer = BatchSerializer(batch)
        return Response(serializer.data)

    def delete(self, request, pk, format=None):
        batch = self.get_object(pk)
        batch.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class ProfileList(generics.ListCreateAPIView):
    queryset = Profile.objects.all()
    serializer_class = ProfileSerializer


class ProfileDetail(APIView):
    def get_object(self, pk):
        try:
            return Profile.objects.get(pk=pk)
        except Profile.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        profile = self.get_object(pk)
        serializer = ProfileSerializer(profile)
        return Response(serializer.data)

    def delete(self, request, pk, format=None):
        profile = self.get_object(pk)
        profile.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class DendrogramList(generics.ListCreateAPIView):
    queryset = Dendrogram.objects.all()
    serializer_class = DendrogramSerializer


class DendrogramDetail(APIView):
    def get_object(self, pk):
        try:
            return Dendrogram.objects.get(pk=pk)
        except Profile.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        dendrogram = self.get_object(pk)
        serializer = DendrogramSerializer(dendrogram)
        return Response(serializer.data)

    def put(self, request, pk, format=None):
        dendrogram = self.get_object(pk)
        serializer = DendrogramSerializer(dendrogram, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk, format=None):
        dendrogram = self.get_object(pk)
        dendrogram.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class Plotting(APIView):
    def post(self, request, format=None):
        serializer = PlotingSerializer(data=request.data)
        if serializer.is_valid():
            plot_dendrogram.delay(str(serializer.data["id"]), serializer.data["linkage"])
            return Response(serializer.data, status=status.HTTP_202_ACCEPTED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
