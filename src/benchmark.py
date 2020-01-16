import collections
# from typing import Dict, Iterable
import click
import ms_io
import logging
import metrics as mx
import pyopenms
import spectrum_utils.spectrum as sus


logger = logging.getLogger('specpride')
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command('evaluate_representatives',
    short_help='Compare an mgf with clustered spectra with a representatives mgf by metric x'
)
@click.option(
    '--cluster_members_file_name', '-c',
    help='MGF file containing the cluster member spectra',
    required=True
)
@click.option(
    '--representatives_file_name', '-r',
    help='MGF file containing the representative spectra',
    required=True
)
@click.option(
    '--out_file_name', '-o',
    help='JSON file name that saves the distances',
    required=True
)
@click.option(
    '--cluster_members_identifications_file_name', '-i',
    help='idXML file containing the identifications for cluster member spectra',
)
@click.option(
    '--representatives_identifications_file_name', '-j',
    help='idXML file containing the identifications for representative spectra',
)
@click.option(
    '--metric_option', '-m',
    help='Distance metric that will be used',
    default="average_cos_dist",
    show_default=True
)
@click.option(
    '--verbose', '-v',
    help='Print verbose logging messages',
    default=False
)
def evaluate_representatives(
    cluster_members_file_name,
    representatives_file_name,
    out_file_name,
    cluster_members_identifications_file_name,
    representatives_identifications_file_name,
    metric_option="average_cos_dist",
    verbose=False,
):
    if metric_option == "average_cos_dist":
        metric = mx.average_cos_dist
    elif metric_option == "fraction_of_by":
        metric = mx.fraction_of_by
        if representatives_identifications_file_name is None:
            raise ValueError("No representatives identifications provided")
        if representatives_identifications_file_name is None:
            raise ValueError("No representatives identifications provided")
    else:
        raise ValueError("Metric not supported")
    logging.info(f'Reading cluster member spectra from {cluster_members_file_name}')
    cluster_member_spectra = ms_io.read_cluster_spectra(
        cluster_members_file_name
    )
    if cluster_members_identifications_file_name is not None:
        cluster_members_identifications = ms_io.read_idXML(
            cluster_members_identifications_file_name
        )
        cluster_member_spectra = annotate(
            cluster_member_spectra,
            cluster_members_identifications
        )
    clusters = collections.defaultdict(list)
    logging.info(f"Grouping cluster member spectra")
    for cluster_member in cluster_member_spectra.values():
        clusters[cluster_member.cluster].append(cluster_member)
    logging.info(f"Reading representatative member spectra from {representatives_file_name}")
    representative_spectra = ms_io.read_cluster_spectra(
        representatives_file_name,
        usi_present=False
    )
    if representatives_identifications_file_name is not None:
        representatives_identifications = ms_io.read_idXML(
            representatives_identifications_file_name
        )
        representative_spectra = annotate(
            representative_spectra,
            representatives_identifications
        )
    logging.info(f"Calculating distances")
    distances = {}
    i = 0
    for cluster, cluster_members in clusters.items():
        # TODO: Quick testing purpose, remove 3 lines for full analysis
        i += 1
        if i > 100:
            break
        if cluster not in representative_spectra:
            if verbose:
                logging.info(f"Cluster: {cluster} has no representative spectrum")
            continue
        representative_spectrum = representative_spectra[cluster]
        distance = metric(representative_spectrum, cluster_members)
        if verbose:
            logging.info(f"Cluster: {cluster}, Distance: {distance}")
        distances[cluster] = distance
    logging.info(f"Saving to JSON file {out_file_name}")
    ms_io.write_distance_dict_to_json(out_file_name, distances)


def annotate(spectra, identifications):
    # TODO: Currently idXML returns index based spectrum reference instead of
    # title based spectrum reference
    new_spectra = {}
    for spectrum_index, (spectrum_id, spectrum) in enumerate(
        sorted(spectra.items()),
        1
    ):
        identification = identifications[spectrum_index]
        sequence = pyopenms.AASequence().fromString(identification)
        modifications = {}
        for residue_index, residue in enumerate(sequence):
            if residue.isModified():
                delta_mass = residue.getModification().getDiffMonoMass()
                modifications[residue_index] = delta_mass
        new_spectrum = sus.spectrum.MsmsSpectrum(
            spectrum.identifier,
            spectrum.precursor_mz,
            spectrum.precursor_charge,
            spectrum.mz,
            spectrum.intensity,
            peptide=sequence.toUnmodifiedString(),
            modifications=modifications
        )
        new_spectra[spectrum_id] = new_spectrum
    return new_spectra


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


cli.add_command(evaluate_representatives)


if __name__ == '__main__':
    try:
        logging.basicConfig(format='{asctime} [{levelname}/{processName}] '
                                   '{module}.{funcName} : {message}',
                            style='{', level=logging.DEBUG, force=True)
    except ValueError:
        # force argument not supported on python 3.6
        logging.basicConfig(format='{asctime} [{levelname}/{processName}] '
                                   '{module}.{funcName} : {message}',
                            style='{', level=logging.DEBUG)

    cli()

    logging.shutdown()
