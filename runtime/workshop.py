"""Small, inspectable research calculations for the three course pilots.

No model output is manufactured here. GeoGPT dialogue is a separate browser step.
"""
from pathlib import Path
import io
import os
import time
import json
from numbers import Real
import numpy as np
import pandas as pd
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parent / '.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve, hilbert, welch
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import train_test_split

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'outputs'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})

def save(fig, name):
    fig.savefig(OUT / (name + '.png'), dpi=160, bbox_inches='tight')
    fig.savefig(OUT / (name + '.svg'), bbox_inches='tight')
    plt.close(fig)

def fossils(rate_m_per_ma=20., event_allowance_ma=0., allow_reworking=False):
    """One-sided occurrence constraints; sensitivity bounds, not probability CIs."""
    started = time.perf_counter()
    if rate_m_per_ma <= 0 or not 0 <= event_allowance_ma <= 1:
        raise ValueError('Positive sedimentation rate and allowance between 0 and 1 Ma required')
    raw = (ROOT / 'data/raw/fossils/M0027A-forams.tsv').read_text()
    df = pd.read_csv(io.StringIO(raw.split('*/', 1)[1].strip()), sep='\t')
    target = df[df['Depth sed [m]'].between(225, 310)].copy()
    # Published GTS2004 event column, IODP 342 methods Table T2.
    # An in situ occurrence cannot predate its species' first occurrence or
    # postdate its last occurrence. Local absences impose NO constraint here.
    younger = target[target['P. sicana'].notna()]['Depth sed [m]'].max()
    older = target[target['C. dissimilis'].notna()]['Depth sed [m]'].min()
    dz = older - younger
    min_elapsed = max(0., 17.54 - 16.38 - 2 * event_allowance_ma)
    if allow_reworking:
        min_elapsed = 0.  # The old specimen no longer supplies a depositional lower age bound.
    min_hiatus = max(0., min_elapsed - dz / rate_m_per_ma)
    rates = np.linspace(2, 60, 150)
    fig, ax = plt.subplots(1, 3, figsize=(12, 4.5), layout='constrained')
    for x, col in enumerate(['P. sicana', 'C. dissimilis']):
        z = target.loc[target[col].notna(), 'Depth sed [m]']
        ax[0].scatter(np.full(len(z), x), z, s=45)
    ax[0].axhspan(younger, older, alpha=.12, color='#cd8c26')
    ax[0].set(xticks=[0,1], xticklabels=['P. sicana','C. dissimilis'], ylabel='Published sample depth (mbsf)', title='Observed occurrences')
    ax[0].invert_yaxis()
    ax[1].plot(rates, np.maximum(0, min_elapsed-dz/rates), color='#176e87')
    ax[1].scatter([rate_m_per_ma], [min_hiatus], color='#b34a30')
    ax[1].set(xlabel='Assumed preserved accumulation (m/Ma)', ylabel='Minimum hiatus required (Ma)', title='Conditional on rate and in situ fossils')
    sr = pd.read_csv(ROOT / 'data/prepared/sr-context.csv')
    ax[2].scatter(sr['age_ma'], sr['depth_mcd'], color='#665a90')
    ax[2].set(xlabel='Published Sr age estimate (Ma)', ylabel='Published depth (mcd)', title='Separate Sr evidence: datum not harmonized')
    ax[2].invert_yaxis()
    fig.suptitle('M0027A: missing time or slow accumulation?')
    save(fig, '01-fossil-constraints')
    target.to_csv(OUT / '01-occurrence-evidence.csv', index=False)
    result = dict(shallower_sample_mbsf=float(younger), deeper_sample_mbsf=float(older), thickness_m=float(dz), minimum_elapsed_ma=float(min_elapsed), assumed_rate_m_per_ma=rate_m_per_ma, minimum_conditional_hiatus_ma=float(min_hiatus), event_allowance_ma=event_allowance_ma, allow_reworking=allow_reworking, conclusion='Hiatus is rate- and taphonomy-dependent; presence data alone do not require one.', sr_depths_harmonized=False, seconds=time.perf_counter()-started)
    pd.DataFrame([result]).to_csv(OUT / '01-results.csv', index=False)
    return result

def load_seismic(time_depth='provided_model'):
    p = ROOT / 'data/raw/f3/data'
    s = (p / 'All_wells_RawData/Lasfiles/F02-1_logs.las').read_text()
    a = np.loadtxt(io.StringIO(s.split('~A',1)[1].split('\n',1)[1]))
    valid = (a[:,1] > 0) & (a[:,2] > 0)
    a = a[valid]
    assert np.all(np.diff(a[:,0]) > 0)
    # LAS declares RHOB kg/m3, DT microseconds/m. Do NOT use a us/ft conversion.
    imp = a[:,1] * 1e6 / a[:,2]
    td = np.loadtxt(p / 'All_wells_RawData/DT_model/F02-1_TD.txt')
    tlog = np.interp(a[:,0], td[:,1], td[:,0]) / 1000
    if time_depth == 'checkshot':
        cs = np.loadtxt(p / 'All_wells_RawData/Checkshot/F02-1_TD.txt')
        tlog = np.interp(a[:,0], cs[:,0], cs[:,1])
    elif time_depth != 'provided_model':
        raise ValueError('Choose provided_model or checkshot')
    section = np.loadtxt(p / 'export_inline362.ascii')
    t = np.arange(section.shape[1]-2) * .004
    well = np.flatnonzero(section[:,1] == 336)
    assert len(well) == 1
    obs = section[well[0],2:]
    ti = np.arange(np.ceil(tlog.min()*1000)/1000, np.floor(tlog.max()*1000)/1000, .001)
    zi = np.interp(ti, tlog, imp)
    r = np.r_[0, np.diff(zi)/(zi[1:]+zi[:-1])]
    return t, obs, ti, r, section, a, tlog

def seismic(max_frequency_hz=55, fixed_phase_deg=None, time_depth='provided_model'):
    """Fit wavelet and a small bulk shift on calibration; freeze for test."""
    started = time.perf_counter()
    if not 15 <= max_frequency_hz <= 80:
        raise ValueError('Frequency cap must be between 15 and 80 Hz')
    t, obs, ti, r, section, logs, tlog = load_seismic(time_depth)
    train = (t >= .55) & (t <= .95)
    test = (t >= 1.10) & (t <= 1.40)
    assert tlog.min() < .55 and tlog.max() > 1.4
    wave_t = np.arange(-.128, .1281, .001)
    phases = np.arange(-180,180,30) if fixed_phase_deg is None else [fixed_phase_deg]
    trials = []
    best = None
    for f in np.arange(15, max_frequency_hz+1, 5):
        q = np.pi*f*wave_t
        w = (1-2*q*q)*np.exp(-q*q)
        h = np.imag(hilbert(w))
        for phase in phases:
            phi = np.deg2rad(phase)
            syn = fftconvolve(r, np.cos(phi)*w - np.sin(phi)*h, mode='same')
            for shift in np.arange(-.020,.0201,.004):
                x = np.interp(t-shift, ti, syn, left=0, right=0)
                design = np.c_[x[train], np.ones(train.sum())]
                scale, offset = np.linalg.lstsq(design, obs[train], rcond=None)[0]
                # Avoid duplicating a phase reversal through negative scale.
                if scale <= 0:
                    continue
                pred = scale*x + offset
                mse = np.mean((obs[train]-pred[train])**2)
                trials.append((f,phase,shift,mse))
                if best is None or mse < best[0]: best=(mse,f,phase,shift,pred)
    _, freq, phase, shift, pred = best
    rows=[]
    baseline = np.mean(obs[train])
    for label, mask in [('calibration',train),('held_out',test)]:
        rmse = np.sqrt(np.mean((obs[mask]-pred[mask])**2))
        ref = np.sqrt(np.mean((obs[mask]-baseline)**2))
        rows.append(dict(interval=label, correlation=float(np.corrcoef(obs[mask],pred[mask])[0,1]), rmse=float(rmse), mean_baseline_rmse=float(ref), rmse_ratio=float(rmse/ref)))
    fig, ax = plt.subplots(1,3,figsize=(13,4.8),layout='constrained')
    crop = (section[:,1] >= 300) & (section[:,1] <= 450)
    vmax = np.percentile(np.abs(section[crop,2:]),98)
    ax[0].imshow(section[crop,2:].T, extent=[300,450,t[-1],0], aspect='auto', cmap='RdBu_r', vmin=-vmax,vmax=vmax)
    ax[0].axvline(336,color='black',lw=1)
    ax[0].set(ylim=(1.5,.4),xlabel='Crossline',ylabel='Two-way time (s)',title='F3 inline 362 — filtered SEG export')
    for aa,(label,mask) in zip(ax[1:], [('Calibration',train),('Held-out interval',test)]):
        aa.plot(obs[mask],t[mask],color='#333333',label='Observed')
        aa.plot(pred[mask],t[mask],color='#c35b39',label='Synthetic')
        aa.set(xlabel='Amplitude (export units)',ylabel='Two-way time (s)',title=label)
        aa.invert_yaxis();aa.legend(fontsize=8)
    fig.suptitle(f'F02-1: measured density × sonic velocity | {freq:g} Hz, {phase:g}°, {shift*1000:g} ms')
    save(fig,'02-well-tie')
    pd.DataFrame(rows).to_csv(OUT/'02-metrics.csv',index=False)
    pd.DataFrame(trials,columns=['frequency_hz','phase_deg','shift_s','calibration_mse']).to_csv(OUT/'02-search.csv',index=False)
    pd.DataFrame(dict(time_s=t,observed=obs,predicted=pred,calibration=train,held_out=test)).to_csv(OUT/'02-traces.csv',index=False)
    return dict(time_depth=time_depth,frequency_hz=float(freq),phase_deg=float(phase),shift_ms=float(shift*1000),metrics=rows,seconds=time.perf_counter()-started,limitation='Bulk wavelet calibration is not independent stratigraphic validation; no original unfiltered volume is included. Checkshot/model depth-datum reconciliation remains an acceptance gate.')

def resolution(frequency_hz=None, phase_deg=0., target_thickness_m=10.):
    """Resolution experiment grounded in real logs and observed bandwidth.

    The wedge is an explicit idealized experiment, not a claimed observed bed.
    No unvalidated well-to-seismic tie is used to identify a stratigraphic boundary.
    """
    started=time.perf_counter()
    for name, value in [('frequency_hz', frequency_hz), ('phase_deg', phase_deg), ('target_thickness_m', target_thickness_m)]:
        if name == 'frequency_hz' and value is None: continue
        if isinstance(value, bool) or not isinstance(value, Real) or not np.isfinite(value):
            raise ValueError(name + ' must be a finite number')
    if not -180 <= phase_deg <= 180:
        raise ValueError('Use phase between -180 and 180 degrees')
    if float(target_thickness_m) != int(target_thickness_m):
        raise ValueError('Use an integer target thickness from 1 to 60 m; the saved curve has 1 m samples')
    t,obs,ti,r,section,logs,tlog=load_seismic()
    window=(t>=.5)&(t<=1.4)
    ff,power=welch(section[(section[:,1]>=300)&(section[:,1]<=450),2:][:,window],fs=250,nperseg=128,axis=1)
    power=power.mean(axis=0);band=(ff>=10)&(ff<=80)
    observed_peak=float(ff[band][np.argmax(power[band])])
    if frequency_hz is None:frequency_hz=observed_peak
    if not 10<=frequency_hz<=80 or not 1<=target_thickness_m<=60:
        raise ValueError('Use 10–80 Hz and 1–60 m')
    # Adjacent measured-log windows chosen before simulation; medians limit spikes.
    a=logs[(logs[:,0]>=1000)&(logs[:,0]<1010)]
    b=logs[(logs[:,0]>=1010)&(logs[:,0]<1020)]
    vp_a=float(np.median(1e6/a[:,2]));vp_b=float(np.median(1e6/b[:,2]))
    rho_a=float(np.median(a[:,1]));rho_b=float(np.median(b[:,1]))
    rc=(vp_b*rho_b-vp_a*rho_a)/(vp_b*rho_b+vp_a*rho_a)
    # Pad the Hilbert-transform window: short periodic boundaries bias rotated 10 Hz wavelets.
    wave_t=np.arange(-2000,2321,dtype=float)*.0005
    def wave(x):
        q=np.pi*frequency_hz*x
        return (1-2*q*q)*np.exp(-q*q)
    w=wave(wave_t);wh=np.imag(hilbert(w))
    rotated=np.cos(np.deg2rad(phase_deg))*w-np.sin(np.deg2rad(phase_deg))*wh
    thickness=np.arange(1,61,dtype=float)
    traces=np.array([rc*(rotated-np.interp(wave_t-2*h/vp_b,wave_t,rotated,left=0,right=0)) for h in thickness])
    peak_separation=np.abs(wave_t[traces.argmax(axis=1)]-wave_t[traces.argmin(axis=1)])
    apparent=peak_separation*vp_b/2
    tuning=np.ptp(traces,axis=1)/(2*abs(rc))
    apparent_target=float(np.interp(target_thickness_m,thickness,apparent))
    fig,ax=plt.subplots(1,3,figsize=(13,4.3),layout='constrained')
    ax[0].plot(ff,power/power.max(),color='#176e87');ax[0].axvline(frequency_hz,color='#c35b39',ls='--')
    ax[0].set(xlim=(0,100),xlabel='Frequency (Hz)',ylabel='Normalized observed power',title='Filtered F3 section spectrum')
    vmax=np.max(np.abs(traces))
    ax[1].imshow(traces.T,origin='upper',extent=[1,60,wave_t[-1]*1000,wave_t[0]*1000],aspect='auto',cmap='RdBu_r',vmin=-vmax,vmax=vmax)
    ax[1].plot(thickness,2*thickness/vp_b*1000,color='black',ls='--',lw=1)
    ax[1].set(ylim=(100,-40),xlabel='True model bed thickness (m)',ylabel='Time relative to bed top (ms)',title='Idealized wedge — measured contrasts')
    ax[2].plot(thickness,apparent,color='#176e87',label='Peak/trough separation estimate')
    ax[2].plot(thickness,thickness,color='gray',ls='--',label='True thickness')
    ax[2].scatter([target_thickness_m],[apparent_target],color='#c35b39')
    ax[2].set(xlabel='True model bed thickness (m)',ylabel='Apparent thickness (m)',title='Thin-bed interpretation bias');ax[2].legend(fontsize=7)
    fig.suptitle(f'What can this bandwidth resolve? {frequency_hz:.1f} Hz, {phase_deg:g}° | synthetic experiment')
    save(fig,'02-resolution')
    pd.DataFrame(dict(thickness_m=thickness,apparent_peak_trough_thickness_m=apparent,normalized_peak_to_peak_amplitude=tuning)).to_csv(OUT/'02-resolution-curves.csv',index=False)
    result=dict(observed_spectral_peak_hz=observed_peak,assumed_wavelet_frequency_hz=float(frequency_hz),phase_deg=phase_deg,synthetic_dt_s=.0005,wavelet_padding_s=1.,vp_background_m_s=vp_a,vp_bed_m_s=vp_b,rho_background_kg_m3=rho_a,rho_bed_kg_m3=rho_b,reflection_coefficient=float(rc),quarter_wavelength_m=float(vp_b/(4*frequency_hz)),target_true_thickness_m=target_thickness_m,apparent_peak_trough_thickness_m=apparent_target,seconds=time.perf_counter()-started,limitation='Observed spectral peak is not a measured wavelet. Wedge geometry is hypothetical. Quarter wavelength is a heuristic, not a detection threshold or a field-bed measurement.')
    pd.DataFrame([result]).to_csv(OUT/'02-resolution-results.csv',index=False)
    return result

def spatial(holdout='east', buffer_km=3., quality_class=None):
    started=time.perf_counter()
    if holdout not in ['east','north'] or not 1 <= buffer_km <= 10:
        raise ValueError('Choose east/north and buffer 1–10 km')
    df=pd.read_csv(ROOT/'data/prepared/gorkha-grid.csv')
    if quality_class is not None:
        if quality_class not in [0,1]: raise ValueError('Image-quality class must be 0, 1 or None')
        df=df[df.image_quality_class==quality_class]
    df=df.reset_index(drop=True)
    coordinate=df.x_m if holdout=='east' else df.y_m
    threshold=float(coordinate.quantile(.75))
    test=coordinate>=threshold
    train=coordinate < threshold-buffer_km*1000
    assert train.sum()>100 and test.sum()>100
    y=df.landslide_present.to_numpy(dtype=int)
    # This is a full cohort of eligible equal-area cells, not sampled pseudoabsences.
    idx=np.arange(len(df));random_train,random_test=train_test_split(idx,test_size=.25,random_state=42,stratify=y)
    rows=[];predictions={}
    for split,tr,te in [('geographic',np.flatnonzero(train),np.flatnonzero(test)),('random',random_train,random_test)]:
        for name,cols in [('terrain',['slope_mean_deg','relief_m']),('terrain + shaking',['slope_mean_deg','relief_m','pga_pctg'])]:
            X=df[cols].to_numpy()
            pipe=make_pipeline(StandardScaler(),LogisticRegression(C=1.,max_iter=500))
            pipe.fit(X[tr],y[tr]);pred=pipe.predict_proba(X[te])[:,1]
            base=np.full(len(te),y[tr].mean())
            row=dict(split=split,model=name,n_train=len(tr),n_test=len(te),test_prevalence=float(y[te].mean()),average_precision=float(average_precision_score(y[te],pred)),roc_auc=float(roc_auc_score(y[te],pred)),brier=float(brier_score_loss(y[te],pred)),baseline_brier=float(brier_score_loss(y[te],base)))
            rows.append(row)
            if split=='geographic': predictions[name]=(te,pred)
    # Paired spatial-block resampling of held-out errors, conditional on these fits.
    te,p1=predictions['terrain + shaking'];_,p0=predictions['terrain']
    block=(df.x_m.to_numpy()[te]//10000).astype(int)*10000+(df.y_m.to_numpy()[te]//10000).astype(int)
    err0=(y[te]-p0)**2;err1=(y[te]-p1)**2
    ids=np.unique(block);rng=np.random.default_rng(42)
    groups=[np.flatnonzero(block==b) for b in ids]
    improvements=[]
    for _ in range(300):
        chosen=np.concatenate([groups[i] for i in rng.integers(0,len(groups),len(groups))])
        improvements.append(float(np.mean(err0[chosen]-err1[chosen])))
    ci=np.quantile(improvements,[.025,.975]).tolist()
    fig,ax=plt.subplots(1,3,figsize=(13,4.4),layout='constrained')
    ax[0].scatter(df.x_m/1000,df.y_m/1000,c=y,s=3,cmap='Greys',vmin=0,vmax=1)
    ax[0].scatter(df.x_m[test]/1000,df.y_m[test]/1000,facecolors='none',edgecolors='#a65436',s=5,lw=.3)
    ax[0].set(title='Mapped occurrence; holdout outlined',xlabel='UTM easting (km)',ylabel='UTM northing (km)',aspect='equal')
    sc=ax[1].scatter(df.x_m.iloc[te]/1000,df.y_m.iloc[te]/1000,c=p1,s=6,vmin=0,vmax=1,cmap='viridis')
    ax[1].set(title='Withheld-region prediction',xlabel='UTM easting (km)',ylabel='UTM northing (km)',aspect='equal');fig.colorbar(sc,ax=ax[1],label='P(mapped source in 500 m cell)')
    metrics=pd.DataFrame(rows)
    for j,model in enumerate(['terrain','terrain + shaking']):
        sub=metrics[metrics.model==model]
        ax[2].bar(np.arange(2)+j*.32,sub.brier,width=.32,label=model)
    ax[2].set(xticks=[.16,1.16],xticklabels=['Geographic','Random'],ylabel='Brier score (lower is better)',title='Transfer changes the assessment');ax[2].legend(fontsize=8)
    save(fig,'03-spatial-validation');metrics.to_csv(OUT/'03-metrics.csv',index=False)
    df.iloc[te].assign(pred_terrain=p0,pred_terrain_shaking=p1).to_csv(OUT/'03-holdout-predictions.csv',index=False)
    return dict(holdout=holdout,buffer_km=buffer_km,quality_class=quality_class,eligible_cells=len(df),heldout_blocks=len(ids),brier_improvement_95pct_block_bootstrap=ci,metrics=rows,seconds=time.perf_counter()-started,limitation='Within-region earthquake-sequence association; interval conditional on fitted models, not prospective hazard uncertainty. Image-quality codes are not ranked without a verified codebook.')
