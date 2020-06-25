function geoLocFormat(lat, lon) {
  latv = Math.abs(lat);
  latvd = Math.trunc(latv);
  latvm = Math.trunc((latv-latvd)*60);
  latf = latvd+'°'+latvm+'"'+((lat>0)?'N':'S');
  lonv = Math.abs(lon);
  lonvd = Math.trunc(lonv);
  lonvm = Math.trunc((lonv-lonvd)*60);
  lonf = lonvd+'°'+lonvm+'"'+((lon>0)?'E':'W');
  return latf+' '+lonf;
}
$('#geoloc5').leafletLocationPicker({
  alwaysOpen: true,
  mapContainer: "#map"
});
$('#save').click(function () {
  this.storageManager = new LocalStorageManager;
  latlon = document.getElementById('geoloc5').value;
  this.storageManager.setLatLon(latlon);
  location.replace("#mainpage");
});
$('#clearloc').click(function () {
  this.storageManager = new LocalStorageManager;
  this.storageManager.clearLatLon();
  document.getElementById('geoloc5').value = "0,0";
});
this.storageManager = new LocalStorageManager;
var latlon = this.storageManager.getLatLon() || "0,0";
document.getElementById('geoloc5').value = latlon;
ll = latlon.split(",");
var geoLat = ll[0];
var geoLon = ll[1];

var live;
function moonload() {
  live = this.setInterval(function(){$("#moon").load("moondata.html")}, 1000);
}
function moonunload() {
  this.clearInterval(live);
}
function el(id) {
  return document.getElementById(id);
}
function lucheck() {
  var lu = el('liveupdate');
  $('#liveupdate').click(function () {
    if (lu.checked) {
      moonload();
    } else {
      moonunload();
    }
  });
}
moonload();
lucheck(); 

var slive;
var planetarium;
function chg(slive) {
  planetarium = S.virtualsky({
    id: 'starmap',
    projection: 'stereo',
    latitude: parseFloat(geoLat),
    longitude: parseFloat(geoLon),
    ecliptic: true,
    meridian: true,
    gridlines_eq: true,
    showstarlabels: true,
    live: slive
  });
}
function slucheck() {
  var lu = el('sliveupdate');
  $('#sliveupdate').click(function () {
    if (lu.checked) {
      chg(true);
    } else {
      chg(false);
    }
  });
}
slucheck();
$(document).on("pageshow", "#skypage", function() {
  chg(false);
});
