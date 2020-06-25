class CMoonPhase {
  constructor(latitude = null, longitude = null, timezone = null, date = null) {
    if (typeof date !== 'Undefined') {
      date = Date.now()/1000;
    } else if (date instanceof DateTime) {
      date = date.getTimestamp();
    }

    this.timestamp = date;

    this.geoLat = latitude;
    this.geoLon = longitude;
    this.timezone = timezone;

    // Astronomical constants. 1980 January 0.0
    let epoch = 2444238.5;

    // Constants defining the Sun's apparent orbit
    let elonge = 278.833540;    // Ecliptic longitude of the Sun at epoch 1980.0
    let elongp = 282.596403;    // Ecliptic longitude of the Sun at perigee
    let eccent = 0.016718;      // Eccentricity of Earth's orbit
    let sunsmax = 1.495985e8;    // Semi-major axis of Earth's orbit, km
    let sunangsiz = 0.533128;    // Sun's angular size, degrees, at semi-major axis distance

    // Elements of the Moon's orbit, epoch 1980.0
    let mmlong = 64.975464;    // Moon's mean longitude at the epoch
    let mmlongp = 349.383063;    // Mean longitude of the perigee at the epoch
    let mlnode = 151.950429;    // Mean longitude of the node at the epoch
    let minc = 5.145396;      // Inclination of the Moon's orbit
    let mecc = 0.054900;      // Eccentricity of the Moon's orbit
    let mangsiz = 0.5181;      // Moon's angular size at distance a from Earth
    let msmax = 384401;      // Semi-major axis of Moon's orbit in km
    let mparallax = 0.9507;    // Parallax at distance a from Earth
    let synmonth = 29.53058868;  // Synodic month (new Moon to new Moon)

    this.synmonth = synmonth;

    // date is coming in as a UNIX timstamp, so convert it to Julian
    date = date / 86400 + 2440587.5;

    // Calculation of the Sun's position
    var Day = date - epoch;                // Date within epoch
    var N = this.fixangle((360 / 365.2422) * Day);    // Mean anomaly of the Sun
    var M = this.fixangle(N + elonge - elongp);    // Convert from perigee co-ordinates to epoch 1980.0
    var Ec = this.kepler(M, eccent);          // Solve equation of Kepler
    var Ec = Math.sqrt((1 + eccent) / (1 - eccent)) * Math.tan(Ec / 2);
    var Ec = 2 * this.rad2deg(Math.atan(Ec));            // True anomaly
    var Lambdasun = this.fixangle(Ec + elongp);    // Sun's geocentric ecliptic longitude

    var F = ((1 + eccent * Math.cos(this.deg2rad(Ec))) / (1 - eccent * eccent)); // Orbital distance factor
    var SunDist = sunsmax / F;              // Distance to Sun in km
    var SunAng = F * sunangsiz;              // Sun's angular size in degrees

    // Calculation of the Moon's position
    var ml = this.fixangle(13.1763966 * Day + mmlong);        // Moon's mean longitude
    var MM = this.fixangle(ml - 0.1114041 * Day - mmlongp);    // Moon's mean anomaly
    var MN = this.fixangle(mlnode - 0.0529539 * Day);        // Moon's ascending node mean longitude
    var Ev = 1.2739 * Math.sin(this.deg2rad(2 * (ml - Lambdasun) - MM));    // Evection
    var Ae = 0.1858 * Math.sin(this.deg2rad(M));                // Annual equation
    var A3 = 0.37 * Math.sin(this.deg2rad(M));                  // Correction term
    var MmP = MM + Ev - Ae - A3;                  // Corrected anomaly
    var mEc = 6.2886 * Math.sin(this.deg2rad(MmP));                // Correction for the equation of the centre
    var A4 = 0.214 * Math.sin(this.deg2rad(2 * MmP));              // Another correction term
    var lP = ml + Ev + mEc - Ae + A4;                // Corrected longitude
    var V = 0.6583 * Math.sin(this.deg2rad(2 * (lP - Lambdasun)));        // Variation
    var lPP = lP + V;                        // True longitude
    var NP = MN - 0.16 * Math.sin(this.deg2rad(M));              // Corrected longitude of the node
    var y = Math.sin(this.deg2rad(lPP - NP)) * Math.cos(this.deg2rad(minc));      // Y inclination coordinate
    var x = Math.cos(this.deg2rad(lPP - NP));                  // X inclination coordinate

    var Lambdamoon = this.rad2deg(Math.atan2(y, x)) + NP;            // Ecliptic longitude
    var BetaM = this.rad2deg(Math.asin(Math.sin(this.deg2rad(lPP - NP)) * Math.sin(this.deg2rad(minc)))); // Ecliptic latitude

    // Right ascension and declination
    var parts = this.radecazalt(Lambdamoon, BetaM);
    var ra = parts[0];
    var dec = parts[1];
    var az = parts[2];
    var alt = parts[3];

    // Calculation of the phase of the Moon
    var MoonAge = lPP - Lambdasun;                  // Age of the Moon in degrees
    var MoonPhase = (1 - Math.cos(this.deg2rad(MoonAge))) / 2;          // Phase of the Moon

    // Distance of moon from the centre of the Earth
    var MoonDist = (msmax * (1 - mecc * mecc)) / (1 + mecc * Math.cos(this.deg2rad(MmP + mEc)));

    var MoonDFrac = MoonDist / msmax;
    var MoonAng = mangsiz / MoonDFrac;                // Moon's angular diameter
    // MoonPar = mparallax / MoonDFrac;              // Moon's parallax

    // Store results
    this.phase = this.fixangle(MoonAge) / 360;          // Phase (0 to 1)
    this.illumination = MoonPhase;                // Illuminated fraction (0 to 1)
    this.age = synmonth * this.phase;              // Age of moon (days)
    this.distance = MoonDist;                  // Distance (kilometres)
    this.diameter = MoonAng;                    // Angular diameter (degrees)
    this.sundistance = SunDist;                  // Distance to Sun (kilometres)
    this.sundiameter = SunAng;                  // Sun's angular diameter (degrees)
    this.longitude = Lambdamoon;                  // Ecliptic longitude
    this.latitude = BetaM;                  // Ecliptic latitude
    this.rightascension = ra;                  // Right ascension
    this.declination = dec;                  // Declination
    this.azimuth = az;                  // Azimuth
    this.altitude = alt;                  // Altitude
  }

  dhms(hours) {
    var h = Math.trunc(hours);
    var m = Math.trunc((hours-h)*60);
    var s = ((hours-h-m/60)*3600).toFixed(2);
    return [h, m, s];
  }
  getUT() {
    var date = new Date();
    var localhour = date.getHours();
    var localmin = date.getMinutes();
    var localsec = date.getSeconds();
    var dst = date.getTimezoneOffset()/60;
    //var UT =localhour+localmin/60+localsec/3600-dst;
    var UT =localhour+localmin/60+localsec/3600+dst;
    return UT;
  }
  gtojd(gmon, gday, gyear) {
    var yd = (gmon<3)?gyear-1:gyear;
    var md = (gmon<3)?gmon+12:gmon;
    var a = Math.trunc(yd/100);
    var b = ((gyear>1582)+((gyear=1582)*(gmon>10))+((gyear=1582)*(gmon=10)*(gday>=15)))?2-a+Math.trunc(a/4):0;
    var c = (yd<0)?Math.trunc((365.25*yd)-0.75):Math.trunc(365.25*yd);
    var d = Math.trunc(30.6001*(md+1));
    var jd = b+c+d+gday+1720994.5;
    return jd;
  }

  /**
  * Right ascension and declination
  *
  * @param float lamb
  * @param float beta
  * @return string
  */
  radecazalt(eclon, eclat) {
    var date = new Date();
    var Gday = date.getUTCDate();
    var Gmonth = date.getUTCMonth()+1;
    var Gyear = date.getUTCFullYear();
    var jd = this.gtojd(Gmonth, Gday, Gyear);
    
    var t = (jd-2451545)/36525;
    var obliquity = 23.439292 - t*(46.815+t*(0.0006-(t*0.00181)))/3600;
    
    var eclonR = this.deg2rad(eclon);
    var eclatR = this.deg2rad(eclat);
    var obliqR = this.deg2rad(obliquity);
    var dec = Math.asin(Math.sin(eclatR)*Math.cos(obliqR)+Math.cos(eclatR)*Math.sin(obliqR)*Math.sin(eclonR));
    var decD = this.rad2deg(dec);
    var y = Math.sin(eclonR)*Math.cos(obliqR)-Math.tan(eclatR)*Math.sin(obliqR);
    var x = Math.cos(eclonR);
    var ra = Math.atan2(y,x);
    var raD = this.rad2deg(ra);
    var raD = this.fixangle(raD);
    var raH = raD*24/360;
    var parts = this.dhms(raH);
    var rah = parts[0];
    var ram = parts[1];
    var ras = parts[2];
    var parts = this.dhms(decD);
    var decd = parts[0];
    var decm = parts[1];
    var decs = parts[2];
    var raf = rah +"h "+ ram +"m "+ ras +"s";
    var decf = (decd<0)?'-('+Math.abs(decd)+'° '+Math.abs(decm)+'" '+Math.abs(decs)+"')":decd +"° "+ decm +'"'+ decs +"'";
    
    var UT = this.getUT();
    var t = (jd-2451545)/36525;
    var t0 = 6.697374558 + 2400.051336 * t + 0.000025862 * t**2;
    var GST = t0-(24*Math.trunc(t0/24)) + UT*1.002737909;
    var GST = GST-(24*Math.trunc(GST/24));
    var offset = this.geoLon*24/360;
    var LST = GST + offset;
    var LST = LST-(24*Math.trunc(LST/24));
    var h1 = LST - raH;
    var H = (h1<0)?h1+24:h1;
    var HR = this.deg2rad(H*360/24);
    var latR = this.deg2rad(this.geoLat);
    
    var sina = Math.sin(dec)*Math.sin(latR)+Math.cos(dec)*Math.cos(latR)*Math.cos(HR);
    var alt = Math.asin(sina);
    var altD = this.rad2deg(alt);
    var y = -Math.cos(dec)*Math.cos(latR)*Math.sin(HR);
    var x = Math.sin(dec)-Math.sin(latR)*sina;
    var az = Math.atan2(y,x);
    var azD = this.rad2deg(az);
    var azD = this.fixangle(azD);
    var parts = this.dhms(altD);
    var altd = parts[0];
    var altm = parts[1];
    var alts = parts[2];
    var parts = this.dhms(azD);
    var azd = parts[0];
    var azm = parts[1];
    var azs = parts[2];
    var azf = azd +"° "+ azm +'" '+ azs +"'";
    var altf = (altd<0)?'-('+Math.abs(altd)+'° '+Math.abs(altm)+'" '+Math.abs(alts)+"')":altd +"°"+ altm +'"'+ alts +"'";
    return [raf, decf, azf, altf];
  }

  /**
  * Fix angle
  *
  * @param float a
  * @return float
  */
  fixangle(a) {
    return (a - 360 * Math.floor(a / 360));
  }

  /**
  * Kepler
  *
  * @param float m
  * @param float ecc
  * @return float
  */
  kepler(m, ecc) {
    // 1E-6
    let epsilon = 0.000001;
    let e = m = this.deg2rad(m);

    do
    {
      var delta = e - ecc * Math.sin(e) - m;
      e -= delta / (1 - ecc * Math.cos(e));
    }
    while (Math.abs(delta) > epsilon);

    return e;
  }

  /**
  * Calculates time  of the mean new Moon for a given base date.
  *   This argument K to this function is the precomputed synodic month index, given by:
  *   K = (year - 1900) * 12.3685
  *   where year is expressed as a year and fractional year.
  *
  * @param int   date
  * @param float k
  * @return float
  */
  meanphase(date, k) {
    // Time in Julian centuries from 1900 January 0.5
    var jt = (date - 2415020.0) / 36525;
    var t2 = jt * jt;
    var t3 = t2 * jt;
    var nt1 = 2415020.75933 + this.synmonth * k + 0.0001178 * t2 - 0.000000155 * t3 + 0.00033 * Math.sin(this.deg2rad(166.56 + 132.87 * jt - 0.009173 * t2));
    return nt1;
  }

  /**
  * Given a K value used to determine the mean phase of the new moon and a
  *   phase selector (0.0, 0.25, 0.5, 0.75), obtain the true, corrected phase time.
  *
  * @param float k
  * @param float phase
  * @return float|null
  */
  truephase(k, phase) {
    var apcor = false;

    k += phase;        // Add phase to new moon time
    var t = k / 1236.85;      // Time in Julian centuries from 1900 January 0.5
    var t2 = t * t;        // Square for frequent use
    var t3 = t2 * t;        // Cube for frequent use
    var pt = 2415020.75933      // Mean time of phase
      + this.synmonth * k
      + 0.0001178 * t2
      - 0.000000155 * t3
      + 0.00033 * Math.sin(this.deg2rad(166.56 + 132.87 * t - 0.009173 * t2));

    var m = 359.2242 + 29.10535608 * k - 0.0000333 * t2 - 0.00000347 * t3;      // Sun's mean anomaly
    var mprime = 306.0253 + 385.81691806 * k + 0.0107306 * t2 + 0.00001236 * t3;  // Moon's mean anomaly
    var f = 21.2964 + 390.67050646 * k - 0.0016528 * t2 - 0.00000239 * t3;      // Moon's argument of latitude

    if (phase < 0.01 || Math.abs(phase - 0.5) < 0.01)
    {
      // Corrections for New and Full Moon
      pt += (0.1734 - 0.000393 * t) * Math.sin(this.deg2rad(m))
        + 0.0021 * Math.sin(this.deg2rad(2 * m))
        - 0.4068 * Math.sin(this.deg2rad(mprime))
        + 0.0161 * Math.sin(this.deg2rad(2 * mprime))
        - 0.0004 * Math.sin(this.deg2rad(3 * mprime))
        + 0.0104 * Math.sin(this.deg2rad(2 * f))
        - 0.0051 * Math.sin(this.deg2rad(m + mprime))
        - 0.0074 * Math.sin(this.deg2rad(m - mprime))
        + 0.0004 * Math.sin(this.deg2rad(2 * f + m))
        - 0.0004 * Math.sin(this.deg2rad(2 * f - m))
        - 0.0006 * Math.sin(this.deg2rad(2 * f + mprime))
        + 0.0010 * Math.sin(this.deg2rad(2 * f - mprime))
        + 0.0005 * Math.sin(this.deg2rad(m + 2 * mprime));

      apcor = true;
    }
    else if (Math.abs(phase - 0.25) < 0.01 || Math.abs(phase - 0.75) < 0.01)
    {
      pt += (0.1721 - 0.0004 * t) * Math.sin(this.deg2rad(m))
        + 0.0021 * Math.sin(this.deg2rad(2 * m))
        - 0.6280 * Math.sin(this.deg2rad(mprime))
        + 0.0089 * Math.sin(this.deg2rad(2 * mprime))
        - 0.0004 * Math.sin(this.deg2rad(3 * mprime))
        + 0.0079 * Math.sin(this.deg2rad(2 * f))
        - 0.0119 * Math.sin(this.deg2rad(m + mprime))
        - 0.0047 * Math.sin(this.deg2rad(m - mprime))
        + 0.0003 * Math.sin(this.deg2rad(2 * f + m))
        - 0.0004 * Math.sin(this.deg2rad(2 * f - m))
        - 0.0006 * Math.sin(this.deg2rad(2 * f + mprime))
        + 0.0021 * Math.sin(this.deg2rad(2 * f - mprime))
        + 0.0003 * Math.sin(this.deg2rad(m + 2 * mprime))
        + 0.0004 * Math.sin(this.deg2rad(m - 2 * mprime))
        - 0.0003 * Math.sin(this.deg2rad(2 * m + mprime));

      // First and last quarter corrections
      if (phase < 0.25)
      {
        pt += 0.0028 - 0.0004 * Math.cos(this.deg2rad(m)) + 0.0003 * Math.cos(this.deg2rad(mprime));
      }
      else
      {
        pt += -0.0028 + 0.0004 * Math.cos(this.deg2rad(m)) - 0.0003 * Math.cos(this.deg2rad(mprime));
        apcor = true;
      }
    }

    return apcor ? pt : null;
  }

  /**
  * Find time of phases of the moon which surround the current date. Five phases are found, starting and
  *   ending with the new moons which bound the current lunation.
  *
  * @return void
  */
  phasehunt() {
    var sdate = this.utc_to_julian(this.timestamp);
    var adate = sdate - 45;
    var ats = this.timestamp - 86400 * 45;
    var uats = new Date(ats);
    var yy = uats.getFullYear();
    var mm = uats.getMonth()+1;

    var k1 = Math.floor((yy + ((mm - 1) * (1 / 12)) - 1900) * 12.3685);
    var adate =  this.meanphase(adate, k1);
    var nt1 = this.meanphase(adate, k1);

    while (true)
    {
      adate += this.synmonth;
      var k2 = k1 + 1;
      var nt2 = this.meanphase( adate, k2 );

      // If nt2 is close to sdate, then mean phase isn't good enough, we have to be more accurate
      if (Math.abs(nt2 - sdate) < 0.75)
      {
        nt2 = this.truephase(k2, 0.0);
      }

      if (nt1 <= sdate && nt2 > sdate)
      {
        break;
      }

      nt1 = nt2;
      k1 = k2;
    }

    // Results in Julian dates
    var dates = [
      this.truephase(k1, 0.0),
      this.truephase(k1, 0.25),
      this.truephase(k1, 0.5),
      this.truephase(k1, 0.75),
      this.truephase(k2, 0.0),
      this.truephase(k2, 0.25),
      this.truephase(k2, 0.5),
      this.truephase(k2, 0.75)
    ];

    this.quarters = [];
    dates.forEach((gjdate, index)=>{
      this.quarters[index] = this.ts_to_human((gjdate - 2440587.5) * 86400);
    });
  }

  ts_to_human(ts) {
    return new Date(ts*1000).toISOString().replace(/T|Z/g," ");
  }
  /**
  * UTC to Julian
  *
  * @param int ts
  * @return float
  */
  utc_to_julian(timestamp) {
    return timestamp / 86400 + 2440587.5;
  }

  /**
  * Get moon phase
  *
  * @return float
  */
  phase() {
    return this.phase;
  }

  /**
  * Get moon properties
  *
  * @param string property_name
  * @return int|float|array|null
  */
  get(property_name)
  {
    return this[property_name] ?? null;
  }

  array_flip(arr) {
    var key, tmp_ar = {};
    for ( key in arr ) {
      if ( arr.hasOwnProperty( key ) ) {
        tmp_ar[arr[key]] = key;
      }
    }
    return tmp_ar;
  }
  /**
  * Get moon phase data
  *
  * @param string name
  * @return float
  */
  get_phase(name) {
    var phases = [
      'new_moon',
      'first_quarter',
      'full_moon', 
      'last_quarter',
      'next_new_moon',
      'next_first_quarter',
      'next_full_moon',
      'next_last_quarter',
    ];

    if (typeof this.quarters !== 'Undefined') {
      this.phasehunt();
    }
    return this.quarters[this.array_flip(phases)[name]] ?? null;
  }

  /**
  * Get current phase name
  *   There are eight phases, evenly split.
  *   A "New Moon" occupies the 1/16th phases either side of phase = 0, and the rest follow from that.
  *
  * @return string
  */
  phase_name() {
    var names = [
      'New Moon',
      'Waxing Crescent',
      'First Quarter',
      'Waxing Gibbous',
      'Full Moon',
      'Waning Gibbous',
      'Third Quarter',
      'Waning Crescent',
      'New Moon',
    ];

    return names[Math.floor((this.phase + 0.0625) * 8)];
  }
  rad2deg(r) {
    return r*180/Math.PI;
  }
  deg2rad(d) {
    return d*Math.PI/180;
  }
}
