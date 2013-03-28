from django.conf.urls.defaults import patterns, include, url
from django.contrib import admin
from django.conf import settings
from django.http.response import HttpResponseRedirect
admin.autodiscover()


urlpatterns = patterns(
        'django.contrib.staticfiles.views', url(r'^site_media/(?P<path>.*)$', 'serve'),
)

urlpatterns += patterns('django.views.generic.simple',
    # Example:
    # (r'^atlasviewer/', include('atlasviewer.foo.urls')),
    url(r'^$', lambda x: HttpResponseRedirect('/admin/')),
    
    url(r'^admin/', include(admin.site.urls)),
    
    # Transcript app
    url(r'^transcript/', include('glasslab.atlasviewer.transcript.urls')),

)
