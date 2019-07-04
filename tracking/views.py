import json
from rest_framework import generics
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from django.http import Http404
from tracking.serializers import ProfileSerializer, TrackedResultsSerializer
from tracking.tasks import track
from tracking.models import Profile, TrackedResults


class ProfileList(generics.ListCreateAPIView):
    def get(self, request, format=None):
        profile = Profile.objects.all()
        serializer = ProfileSerializer(profile)
        return Response(serializer.data)

    def post(self, request, format=None):
        serializer = ProfileSerializer(data=request.data)
        if serializer.is_valid():
            serializer.save()
            id = str(serializer.data["id"])
            profile = ProfileDetail.get_object(pk=id)
            track.delay(id, str(profile.profile_db))
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


class ProfileDetail(APIView):
    @classmethod
    def get_object(cls, pk):
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


class TrackedResultsList(generics.ListCreateAPIView):
    queryset = TrackedResults.objects.all()
    serializer_class = TrackedResultsSerializer


class TrackedResultsDetail(APIView):
    @classmethod
    def get_object(cls, pk):
        try:
            return TrackedResults.objects.get(pk=pk)
        except TrackedResults.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        results = self.get_object(pk)
        serializer = TrackedResultsSerializer(results)
        data = serializer.data.copy()
        data["json"] = json.loads(results.json.read().decode("utf-8"))
        return Response(data)

    def put(self, request, pk, format=None):
        results = self.get_object(pk)
        serializer = TrackedResultsSerializer(results, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk, format=None):
        results = self.get_object(pk)
        results.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)

