import marimo

__generated_with = "0.23.8"
app = marimo.App(width="medium", auto_download=["ipynb"])


@app.cell
def _():
    import marimo as mo
    import sys
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.time import Time, TimeDelta
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from astropy.visualization import time_support
    from matplotlib.patches import Ellipse

    # 1. Standard environment path setup
    site_pkgs_path = '/opt/anaconda3/envs/radio/lib/python3.14/site-packages'
    if site_pkgs_path not in sys.path:
        sys.path.insert(0, site_pkgs_path)

    import vlbiplanobs.stations as stations
    import vlbiplanobs.sources as sources
    import vlbiplanobs.observation as observation
    import vlbiplanobs.scheduler as scheduler

    return (
        Ellipse,
        SkyCoord,
        Time,
        TimeDelta,
        mdates,
        mo,
        np,
        observation,
        plt,
        scheduler,
        sources,
        stations,
        u,
    )


@app.function
def SED(lamb, lamb0, f0, alpha):
    return f0 * (lamb0/lamb)**alpha


@app.cell
def _(stations):
    master_catalog = stations.Stations()
    vla_codes = master_catalog.filter_networks('VLBA').stations[-2:]
    return (master_catalog,)


@app.cell
def _(mo):
    # Target Information Stuff
    target_name = mo.ui.text(
        placeholder='Search...',
        label='Target',
        value='M87'
    )

    date = mo.ui.date(label="Date")
    clockstart = mo.ui.number(start=0, stop=23.99, value=12.00, step=0.5, label='Time')

    duration = mo.ui.number(start=1, stop=24, step=0.25, label='Duration')

    efficiency = mo.ui.number(start=0, stop=100, step=1, value=70, label='TOS')

    target_flux = mo.ui.number(
        start=0.001,
        stop = 10000,
        label='Flux (mJy)',
        value=1
    )

    target_radius = mo.ui.number(
        start=0.001,
        stop = 10*3600,
        label='Radius (mas)',
        value=1
    )
    return (
        clockstart,
        date,
        duration,
        efficiency,
        target_flux,
        target_name,
        target_radius,
    )


@app.cell
def _(
    SkyCoord,
    Time,
    TimeDelta,
    clockstart,
    date,
    duration,
    np,
    sources,
    target_name,
    u,
):
    target_coords = SkyCoord.from_name(target_name.value)
    target_source = sources.Source(name=target_name.value, coordinates=target_coords)

    start_time = Time(date.value.isoformat() + "T"+str(int(clockstart.value))+":"+ str(int(60*(clockstart.value - int(clockstart.value)))) + ":00", format='isot', scale='utc')

    # 4. Formulate the Scans structure required by the Observation type validation
    obs_time = duration.value * u.hour # Match your exact website slider configuration
    time_steps = np.arange(0, obs_time.to(u.min).value + 15, 1) * u.min
    observation_time_series = start_time + TimeDelta(time_steps)



    single_scan = sources.Scan(source=target_source, duration=obs_time)
    scan_block = sources.ScanBlock(scans=[single_scan])
    scans_dict = {'target': scan_block}
    return observation_time_series, scans_dict


@app.cell
def _(master_catalog, u):
    # Network Stuff

    all_nets = list(master_catalog.get_networks_from_configfile().keys())
    dataratelist = {'4Mbps':2**2 * u.Mbit/u.s,
                    '8Mbps':2**3 * u.Mbit/u.s,
                    '16Mbps':2**4 * u.Mbit/u.s,
                    '32Mbps':2**5 * u.Mbit/u.s,
                    '64Mbps':2**6 * u.Mbit/u.s,
                    '128Mbps':2**7 * u.Mbit/u.s,
                    '256Mbps':2**8 * u.Mbit/u.s,
                    '512Gbps':2**9 * u.Mbit/u.s,
                    '1Gbps':2**10 * u.Mbit/u.s,
                    '2Gbps':2**11 * u.Mbit/u.s,
                    '4Gbps':2**12 * u.Mbit/u.s,
                    '8Gbps':2**13 * u.Mbit/u.s,
                    '16Gbps':2**14 * u.Mbit/u.s,
                    '33Gbps':2**15 * u.Mbit/u.s,
                   }
    subbandlist = list(range(1,9,1))
    channellist = [2**i for i in range(5,15)]
    polarizationlist = {'Full':4, 'Dual':2, 'Single':1}
    return all_nets, channellist, dataratelist, polarizationlist, subbandlist


@app.cell
def _(all_nets, channellist, dataratelist, mo, polarizationlist, subbandlist):
    network = mo.ui.dropdown(
        options=all_nets,
        label='Network',
        value=all_nets[0]
    )




    datarate_opts=mo.ui.dropdown(
        options=dataratelist,
        label='Date Rate',
        value='4Gbps'
    )

    subbands_opts = mo.ui.dropdown(
        options=subbandlist,
        label='Subbands',
        value=subbandlist[2]
    )

    channels_opts = mo.ui.dropdown(
        options=channellist,
        label='Channels',
        value=channellist[3]
    )

    pol_opts = mo.ui.dropdown(
        options=polarizationlist,
        label='Polarization',
        value='Full'
    )
    return channels_opts, datarate_opts, network, pol_opts, subbands_opts


@app.cell
def _(master_catalog, network):
    network_object_master = master_catalog.filter_networks(network.value)
    station_list = network_object_master.stations
    return network_object_master, station_list


@app.cell
def _(mo, network_object_master):
    band_counts = {}

    for station in network_object_master:
        # Iterate through the bands this specific telescope has SEFD data for
        # (Using .sefd.keys() to see the supported bands like '6cm', '1.3cm', etc.)
        for band in station.bands:
            band_counts[band] = band_counts.get(band, 0) + 1

    # 4. Filter out bands supported by at least 2 telescopes (required for a baseline)
    available_bands = [band for band, count in band_counts.items() if count >= 2]

    target_band = mo.ui.dropdown(
        options=available_bands,
        label='Band',
        value=available_bands[0]
    )
    return available_bands, target_band


@app.cell
def _(master_catalog, network):
    network_object = master_catalog.filter_networks(network.value)
    return (network_object,)


@app.cell
def _(mo, station_list):
    selected_stations = mo.ui.multiselect(
        options=station_list,
        label='Stations',
        value=station_list
    )

    return (selected_stations,)


@app.cell
def _(
    channels_opts,
    datarate_opts,
    efficiency,
    network_object,
    observation,
    observation_time_series,
    pol_opts,
    scans_dict,
    selected_stations,
    station_list,
    subbands_opts,
    target_band,
):
    for stats_r in station_list:
        network_object.remove_station(stats_r)

    for stats_a in selected_stations.value:
        network_object.add_station(stats_a)

    # 6. Initialize the core Observation instance
    obs = observation.Observation(
        band=target_band.value,
        stations=network_object,
        scans=scans_dict,
        times=observation_time_series, # Multi-point tracking
        datarate=datarate_opts.value,
        subbands=subbands_opts.value,
        channels=channels_opts.value,
        polarizations=pol_opts.value,               
        ontarget=efficiency.value/100                   
    )
    return (obs,)


@app.cell
def _(mo, np, obs, scheduler, target_flux, target_name, target_radius, u):
    try:
        planned_schedule = scheduler.ObservationScheduler(obs)

        noise_output = planned_schedule.obs.thermal_noise()

        beamshape = planned_schedule.obs.synthesized_beam()

        '''
        if isinstance(noise_output, dict) and target_name.value in noise_output:
            noise = noise_output[target_name.value]
            if noise is not None:
                rms = noise.to(u.uJy/u.beam)
            else:
                rms = 1


        if isinstance(beamshape, dict) and target_name.value in beamshape:
            beam = beamshape[target_name.value]
            if beam is not None:
                bmaj  = beam['bmaj']
                bmin  = beam['bmin']
                PA    = beam['pa']
                barea = bmaj * bmin
            else:
                barea = 1
        '''

        rms = noise_output[target_name.value].to(u.mJy/u.beam)
        beam = beamshape[target_name.value]
        bmaj  = beam['bmaj']
        bmin  = beam['bmin']
        PA    = beam['pa']
        barea = bmaj * bmin


        Nbeams = (target_radius.value * u.mas)**2 / barea

        SNR = ((target_flux.value * u.mJy)/rms / np.sqrt(Nbeams)).value
    except:
        rms = None
        Nbeams = None
        SNR = None
        mo.callout(f"{target_name.value} not observable with this configuration.")
    return Nbeams, PA, SNR, bmaj, bmin, rms


@app.cell
def _(mdates, np, obs, plt, target_name):
    # 3. Create the plot workspace
    fig, ax = plt.subplots(figsize=(10, 6))

    # 4. Loop through your active stations and calculate the elevation tracks
    relative_hours = obs.times.mjd
    elevations = obs.elevations()[target_name.value]
    for stat_e in obs.stations.station_codenames:
        ax.plot(relative_hours, elevations[stat_e], label=stat_e, lw=2)

    # 5. Add formatting matching a professional visibility dashboard
    plt.title(f"Source Elevation Track for {target_name.value}", fontsize=14, fontweight='bold')
    ax.set_xlabel("Time (UTC)", fontsize=12)
    ax.set_ylabel("Elevation Angle (Degrees)", fontsize=12)

    #ax.axhspan(0, 10, color='grey', alpha=0.75) # Typical VLBI limit
    #ax.axhspan(10, 20, color='grey', alpha=0.5) # Typical VLBI limit
    #ax.axhspan(20, 30, color='grey', alpha=0.25) # Typical VLBI limit

    # 2. Define the gradient's Y-boundaries and Colormap
    y_min, y_max = 0,20      # Y-range for the span
    cmap = 'gray'           # Choose any valid Matplotlib colormap

    # 3. Create the smooth gradient background
    # We pass an array of normalized values to imshow and stretch it over the axes
    gradient = np.linspace(0, 1, 256).reshape(-1, 1)
    ax.imshow(gradient, 
              cmap=cmap, 
              aspect='auto', 
              origin='lower', 
              extent=[ax.get_xlim()[0], ax.get_xlim()[1], y_min, y_max], 
              interpolation='bicubic', 
              zorder=0)

    ax.set_ylim(0, 90)
    ax.set_xlim(relative_hours[0], relative_hours[-1])
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    plt.tight_layout()
    return (fig,)


@app.cell
def _(
    available_bands,
    channels_opts,
    datarate_opts,
    efficiency,
    network_object,
    observation,
    observation_time_series,
    pol_opts,
    scans_dict,
    scheduler,
    subbands_opts,
    target_name,
    u,
):
    rms_array = []
    barea_array = []

    for band_check in available_bands:
        obs_hold = observation.Observation(
            band=band_check,
            stations=network_object,
            scans=scans_dict,
            times=observation_time_series, # Multi-point tracking
            datarate=datarate_opts.value,
            subbands=subbands_opts.value,
            channels=channels_opts.value,
            polarizations=pol_opts.value,   
            ontarget=efficiency.value/100                   # Match the website's overhead overhead-efficiency slider
        )
        plan_hold = scheduler.ObservationScheduler(obs_hold)
        noise_output_hold = plan_hold.obs.thermal_noise()
        beamshape_hold = plan_hold.obs.synthesized_beam()
        rms_array.append(noise_output_hold[target_name.value].to(u.mJy/u.beam))
        barea_array.append(beamshape_hold[target_name.value]['bmaj'].to(u.mas)\
                           * beamshape_hold[target_name.value]['bmin'].to(u.mas)
                          )
    rms_array = u.Quantity(rms_array)
    barea_array = u.Quantity(barea_array)
    return barea_array, rms_array


@app.cell
def _(mo):
    S0val = mo.ui.number(
        start=0.001,
        stop=10000,
        value=1,
        label=r"$F_0$ (mJy)"
    )

    lamb0val = mo.ui.number(
        start=0.087,
        stop=92,
        value=2.3,
        label=r'$\lambda_0$ (cm)'
    )

    alphaval = mo.ui.number(
        start=-10,
        stop=10,
        value=0.6,
        label=r'$\alpha$'
    )
    return S0val, alphaval, lamb0val


@app.cell
def _(
    S0val,
    alphaval,
    available_bands,
    barea_array,
    lamb0val,
    np,
    rms_array,
    target_radius,
    u,
):
    waves = u.Quantity([u.Quantity(a) for a in available_bands])
    fluxes = SED(waves,lamb0val.value * u.cm, S0val.value*u.mJy, alphaval.value)
    beamSNR = fluxes /rms_array / np.sqrt((target_radius.value * u.mas)**2 / barea_array)
    return beamSNR, fluxes, waves


@app.cell
def _(SNR, beamSNR, fluxes, obs, plt, u, waves):
    figSNR, (axSNR, axSED) = plt.subplots(2, 1, sharex=True, figsize=(5,7.5))

    axSNR.axhline(SNR, c='r', alpha=0.5)
    axSNR.axvline(obs.wavelength.value, c='r', alpha=0.5)
 

    axSNR.plot(waves, beamSNR, c="k", marker='.')

    to_GHz = lambda x: (x * u.cm).to(u.GHz, equivalencies=u.spectral()).value

    secax = axSNR.secondary_xaxis('top', functions=(to_GHz,to_GHz))

    axSNR.grid(
        color='gray',     
        linestyle='--',    
        linewidth=0.5,     
        alpha=0.7
    )


    secax.set_xlabel("Frequency (GHz)")

    axSNR.set_ylabel('Image SNR')

    axSNR.set_xscale('log')

    axSED.loglog(waves, fluxes, c="k")

    axSED.set_xlabel('Wavelength (cm)')

    axSED.set_ylabel('Flux (mJy)')

    axSED.grid(
        color='gray',
        which='both',
        linestyle='--',    
        linewidth=0.5,     
        alpha=0.7
    )

    plt.subplots_adjust(wspace=0, hspace=0)

    return (figSNR,)


@app.cell
def _(Ellipse, PA, bmaj, bmin, plt, target_radius):
    figbeam, axbeam = plt.subplots(figsize=(5,5))

    beam_ellipse = Ellipse(xy=(0,0), width=bmaj.value*2, height=bmin.value*2, angle=PA.value + 90, 
                      edgecolor='blue', facecolor='none',)

    source_ellipse = Ellipse(xy=(0,0), width=target_radius.value*2, height=target_radius.value*2, angle=0, 
                      edgecolor='black', facecolor='none',)

    axbeam.add_patch(beam_ellipse)
    axbeam.add_patch(source_ellipse)

    axbeam.annotate(text=rf'Beam Dimension: {bmaj.value:0.2}$\times${bmin.value:0.2} mas$^2$', xy=(-1.35*target_radius.value,1.35*target_radius.value))

    axbeam.set_xlim(-1.5*target_radius.value,1.5*target_radius.value)
    axbeam.set_ylim(-1.5*target_radius.value,1.5*target_radius.value)
    axbeam.set_aspect('equal')

    axbeam.legend(['Beam', 'Target'], loc='best')

    plt.tight_layout()
    return (figbeam,)


@app.cell
def _(
    channels_opts,
    clockstart,
    datarate_opts,
    date,
    duration,
    efficiency,
    fig,
    mo,
    network,
    pol_opts,
    selected_stations,
    subbands_opts,
    target_band,
    target_flux,
    target_name,
    target_radius,
):
    output = mo.hstack([mo.vstack([target_name,
                                   network,
                                   selected_stations,
                                   target_band,
                                   date,
                                   clockstart,
                                   duration,
                                   efficiency,
                                   target_flux,
                                   target_radius,
                                   datarate_opts,
                                   subbands_opts,
                                   channels_opts,
                                   pol_opts,   
                                ]),
                            fig
                           ])
    return (output,)


@app.cell
def _(output):
    output
    return


@app.cell
def _(Nbeams, S0val, SNR, alphaval, figSNR, figbeam, lamb0val, mo, rms, u):
    mo.hstack([
        mo.vstack([
            mo.md("""
                    ## For the current observation configuration:
            """),
            mo.md(rf"""
                    The predicted thermal noise is $\sigma\approx$ {rms.to(u.uJy/u.beam):0.2f}
            """),
            mo.md(rf"""
                    There will be $N\approx$ {Nbeams:0.0f} Beams across the objects
            """),
            mo.md(rf"""
                    The image will have SNR $\approx$ {SNR:0.2f}
            """),
            figbeam
        ]),
        mo.vstack([
            mo.md("""
                    ## User input SED Information
            """),
            S0val,
            lamb0val,
            alphaval,
            mo.mpl.interactive(figSNR)
        ])
    ])
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
