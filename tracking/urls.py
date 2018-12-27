from django.urls import path

from . import views

urlpatterns = [
    path('sequence/', views.SequenceList.as_view(), name="sequence-list"),
    path('sequence/<uuid:pk>/', views.SequenceDetail.as_view(), name="sequence-detail"),
    path('results/', views.TrackedResultsList.as_view(), name="results-list"),
    path('results/<uuid:pk>/', views.TrackedResultsDetail.as_view(), name="results-detail"),
    path('tracking/', views.Tracking.as_view(), name="tracking"),
]
