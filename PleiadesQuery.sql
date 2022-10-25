select
    gs.source_id, gs.ra, gs.ra_error, gs.dec, gs.dec_error, gs.parallax, gs.parallax_error, gs.pmra, gs.pmra_error, gs.pmdec, gs.pmdec_error,
    gs.radial_velocity, gs.radial_velocity_error,
    gs.phot_g_mean_mag, gs.phot_bp_mean_mag, gs.phot_rp_mean_mag,
    ap.mh_gspspec, ap.mh_gspspec_lower, ap.mh_gspspec_upper, ap.fem_gspspec, ap.fem_gspspec_lower, ap.fem_gspspec_upper,
    ap.ag_gspphot, ap.ag_gspphot_lower, ap.ag_gspphot_upper, ap.mg_gspphot, ap.mg_gspphot_lower, ap.mg_gspphot_upper,
    cep.PF, cep.PF_ERROR, cep.INT_AVERAGE_G, cep.INT_AVERAGE_G_ERROR
from
    gaiadr3.gaia_source as gs
join gaiadr3.astrophysical_parameters as ap using (source_id)
left join gaiadr3.vari_cepheid as cep using (source_id)
where
    (gs.phot_g_mean_mag is not Null)
    and (gs.phot_bp_mean_mag is not Null)
    and (gs.phot_rp_mean_mag is not Null)
    and (1 = CONTAINS(
        POINT(56.85, 24.1167),
        CIRCLE(gs.ra, gs.dec, 120./60.)))
