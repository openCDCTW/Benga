from django.urls import path
from . import views

urlpatterns = [
    path('', views.IndexView.as_view()),
    path('non-release/', views.NonRelease.as_view()),
]
