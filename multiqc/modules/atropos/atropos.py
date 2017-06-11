#!/usr/bin/env python

""" MultiQC module to parse output from atropos. """

from __future__ import print_function
import logging
import re
from distutils.version import StrictVersion
from collections import OrderedDict

from multiqc import config
from multiqc.plots import linegraph, bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    Atropos module class, parses text logs. Adapted from atropos class.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Atropos', anchor='atropos',
        href='https://atropos.readthedocs.io/',
        info="is a tool to find and remove adapter sequences, primers, poly-A"\
         "tails and other types of unwanted sequence from your high-throughput"\
         " sequencing reads.")

        # Find and load any Atropos reports
        self.atropos_data = dict()
        self.atropos_length_counts = dict()
        self.atropos_length_exp = dict()
        self.atropos_length_obsexp = dict()
        self.atropos_found_adapters = []

        for f in self.find_log_files('atropos', filehandles=True):
            self.parse_atropos_logs(f)

        # Filter to strip out ignored sample names
        self.atropos_data = self.ignore_samples(self.atropos_data)

        if len(self.atropos_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.atropos_data)))

        # Write parsed report data to a file
        self.write_data_file(self.atropos_data, 'multiqc_atropos')

        # Basic Stats Table
        self.atropos_general_stats_table()

        # Trimming Length Profiles
        self.atropos_length_trimmed_plot()

        # Bar plot
        self.atropos_barplot()


    def parse_atropos_logs(self, f):
        """ Go through log file looking for atropos output """
        fh = f['f']

        # adapted from cutadapt version adjustment. 
        # if newer version released, can update similarly.
        regexes = {
            '1.1.5_paired': {
                'bp_processed': "Total bp processed:\s*([\d,]+)\s",
                'bp_written': "Total bp written \(filtered\):\s*([\d,]+)\s",
                'quality_trimmed': "Quality-trimmed:\s*([\d,]+)\s",
                'r_processed': "Total read pairs processed:\s*([\d,]+)\s",
                'r1_with_adapters': "Read 1 with adapter:\s*([\d,]+)\s",
                'r2_with_adapters': "Read 2 with adapter:\s*([\d,]+)\s",
                'r_written': "Pairs written \(passing filters\):\s*([\d,]+)\s"
            },
            '1.1.5_single': {
                'bp_processed': "Total bp processed:\s*([\d,]+)\s",
                'bp_written': "Total bp written \(filtered\):\s*([\d,]+)\s",
                'quality_trimmed': "Quality-trimmed:\s*([\d,]+)\s",
                'r_processed': "Total reads processed:\s*([\d,]+)\s",
                'r_with_adapters': "Reads with adapter:\s*([\d,]+)\s",
                'r_written': "Reads written \(passing filters\):\s*([\d,]+)\s"
            }
        }
        s_name = None
        atropos_version = '1.1.5'
        paired = 'single'

        for l in fh:
            # New log starting
            if 'Atropos version' in l:
                s_name = None
                a_version = re.match(r'Atropos version: ([\d\.]+)', l)

                # This code does nothing now, but can update to allow for 
                # alternative revex versioning
                if a_version:
                    try:
                        assert(StrictVersion(a_version.group(1)) <=
                                             StrictVersion('1.1.5'))
                        atropos_version = '1.1.5'
                    except:
                        atropos_version = '1.1.5'

            # Get mode (single or paired)
            if l.startswith('Input format: '):
                if 'Paired' in l:
                    paired = 'paired'

            # Get sample name
            if l.startswith('Sample ID: '):
                s_name = l.split()[-1]
                # Manage case where sample name is '-' (reading from stdin)
                if s_name == '-':
                    s_name = f['s_name']
                else:
                    s_name = self.clean_s_name(s_name, f['root'])
                if s_name in self.atropos_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.atropos_data[s_name] = dict()

            if s_name is not None:
                self.add_data_source(f, s_name)

                a_regex_type = '%s_%s' % (atropos_version, paired)

                # Search regexes for overview stats
                for k, r in regexes[a_regex_type].items():
                    match = re.search(r, l)
                    if match:
                        self.atropos_data[s_name][k] = int(match.group(1).replace(',', ''))

                # Histogram showing lengths trimmed for each adapter
                if re.match( r'.*Adapter \d', l):

                    adapter_num = l.rstrip()

                    next(fh) # skip the ----
                    next(fh) # skip the blank line
                    next(fh) # header
                    next(fh) # skip the ----

                    l = next(fh)

                    adapter_seq = l.split()[0]

                    if adapter_seq not in self.atropos_found_adapters:
                        self.atropos_found_adapters.append(adapter_seq)

                    # initialize the hist dict for that adapter
                    if adapter_seq not in self.atropos_length_counts:
                        self.atropos_length_counts[adapter_seq] = dict()
                        self.atropos_length_exp[adapter_seq] = dict()

                    self.atropos_length_counts[adapter_seq][s_name] = dict()
                    self.atropos_length_exp[adapter_seq][s_name] = dict()

                    at_header = False
                    in_hist = False
                    
                    # Nested loop to read this section while the regex matches
                    for l in fh:
                        # if you're in the histogram portion, read it
                        if in_hist:
                            r_seqs = re.search("^\s+(\d+)\s+([\d,]+)\s+([\d\.,]+)", l)
                            if r_seqs:
                                a_len = int(r_seqs.group(1))
                                self.atropos_length_counts[adapter_seq][s_name][a_len] = int(r_seqs.group(2).replace(',', ''))
                                self.atropos_length_exp[adapter_seq][s_name][a_len] = float(r_seqs.group(3).replace(',', ''))
                            else:
                                break
                        
                        # otherwise see if you've already hit the header
                        elif at_header:
                            # if you've seen the header and are at the divider
                            if l.startswith('---'):
                                # you've hit the numbers
                                in_hist = True
                                continue
                            else:
                                # otherwise just wait
                                continue

                        # otherwise wait until you get to the histogram header
                        elif l.startswith('Overview of removed sequences:'):
                            at_header = True
                            continue

        # Calculate a few extra numbers of our own
        for s_name, d in self.atropos_data.items():
            self.atropos_data[s_name]['r_lost'] = int(d['r_processed']) - int(d['r_written'])

            if 'bp_processed' in d and 'bp_written' in d:
                self.atropos_data[s_name]['percent_trimmed'] = (float(d['bp_processed'] - d['bp_written']) / d['bp_processed']) * 100
            elif 'bp_processed' in d and 'bp_trimmed' in d:
                self.atropos_data[s_name]['percent_trimmed'] = ((float(d.get('bp_trimmed', 0)) + float(d.get('quality_trimmed', 0))) / d['bp_processed']) * 100


    def atropos_general_stats_table(self):
        """ Take the parsed stats from the atropos report and add it to the
        basic stats table at the top of the report """

        headers = {}
        headers['percent_trimmed'] = {
            'title': '% Trimmed',
            'description': '% Total Base Pairs trimmed',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlBu-rev'
        }
        self.general_stats_addcols(self.atropos_data, headers)


    def atropos_length_trimmed_plot (self):
        """ Generate the trimming length plot """

        for i, adapter in enumerate(self.atropos_found_adapters):
            description = 'This plot shows the number of sequences trimmed by each \
                           length for a given adapter'
            # Adapters found: \n%s' % '\n'.join(['%s: %s' % (x + 1,y) for (x,y) in
            #                                   zip(range(len(self.atropos_found_adapters)),
            #                                       self.atropos_found_adapters)])

            pconfig = {
                'id': 'atropos_plot_length_%s' % i,
                'title': 'Lengths of Trimmed Sequences for adapter: %s' % adapter,
                'ylab': 'Counts',
                'xlab': 'Length Trimmed (bp)',
                'xDecimals': False,
                'ymin': 0,
                'tt_label': '<b>{point.x} bp trimmed</b>: {point.y:.0f}',
                'data_lables': [{'name': adapter, 'ylab': 'Count'}]
            }

            self.add_section(
                name = 'Trimmed length: %s' % adapter,
                anchor = adapter,
                description = description,
                plot = linegraph.plot([self.atropos_length_counts[adapter]], dict(pconfig))
            )



    def atropos_barplot (self):
        """ Make the HighCharts HTML to plot the reads remaining """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['r_written'] =              { 'color': '#437bb1', 'name': 'Surviving Reads' }
        keys['r_lost'] =                { 'color': '#7f0000', 'name': 'Dropped' }

        # Config for the plot
        pconfig = {
            'id': 'atropos_bar_plot',
            'title': 'Atropos',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        self.add_section(name = 'Reads remaining per sample',
                         anchor = 'atropos_reads_remaining_barplot',
                         plot = bargraph.plot(self.atropos_data, keys, pconfig) )

