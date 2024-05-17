
from string import Template
import os

HTML_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'template')

class Report():


    """
    Report object
    """

    def __init__(self, report=''):
        """
        Parameters
        ----------
        report : str
            Report in html
        """
        self.report = report


    def write_params(self, params):
        """
        Parameters
        ----------
        params : Namespace
            Parameters
        """

        with open(os.path.join(HTML_DIR, 'params.html'), 'r') as f:
            template = Template(f.read())
        
        params_report = template.substitute(threshold_c=params.threshold_c,
                                            threshold_r=params.threshold_r,
                                            num_partitions=params.num_partitions,
                                            redundancy_reduction=bool(params.threshold_r),
                                            cluster_pruning=params.prune)
        
        self.report += params_report

    # def write_result(self, results, zip_file):
    #     """
    #     Parameters
    #     ----------
    #     result : list
    #         List of list of clustering results under different thresholds
    #     """

    #     self.report += "<h2>Results</h2>\n"

    #     html_table = '<table>\n'
        
    #     # Add a row to the HTML table for each row in the CSV file
    #     for i, row in enumerate(results):
    #         if i == 0:  # Assuming the first row is the header
    #             html_table += '  <thead>\n    <tr>' + ''.join([f'<th>{cell}</th>' for cell in row]) + '</tr>\n  </thead>\n  <tbody>\n'
    #         else:
    #             html_table += '    <tr>' + ''.join([f'<td>{cell}</td>' for cell in row]) + '</tr>\n'
    #     html_table += '  </tbody>\n</table>\n'

    #     self.report += html_table

    #     # write the donwload button to link
    #     self.report += f"<a href=\"{zip_file}\" download>Download zip results</a>\n"
    
    def write_results(self, results, zip_file):
        """
        Parameters
        ----------
        results : list
            List of clustering results
        zip_file : str
            Path to the zip file
        """

        html_table = '<table>\n'
        have_na = False
        
        # Add a row to the HTML table for each row in the CSV file
        for i, row in enumerate(results):
            if i == 0:  # Assuming the first row is the header
                html_table += '  <thead>\n    <tr>' + ''.join([f'<th>{cell}</th>' for cell in row]) + '</tr>\n  </thead>\n  <tbody>\n'
            elif row[-1] == 'NA':
                html_table += '    <tr>' + ''.join([f'<td>{cell}</td>' for cell in row]) + f'</tr>\n'
                have_na = True
            else:
                html_table += '    <tr>' + ''.join([f'<td>{cell}</td>' for cell in row[:-1]]) + f'<td><a href="{row[-1]}" download>Link</a></td></tr>\n'
        html_table += '  </tbody>\n</table>\n'

        if have_na:
            html_table += "<br>\n* NA indicates the size of maximum cluster exceeds the maximum partition capacity.\n<br>\n"

        with open(os.path.join(HTML_DIR, 'results.html'), 'r') as f:
            template = Template(f.read())
        
        params_report = template.substitute(results_table=html_table,
                                            zip_file=zip_file)
        
        self.report += params_report



    def write_figures(self, figures):
        """
        Parameters
        ----------
        figures : list
            List of figure paths: [(threshold, fig_path_1, fig_path_2, output_file*, sizebar_file*), ...]
        """

        self.report += "<h2>Analysis</h2>\n<hr>\n"

        # load the template
        with open(os.path.join(HTML_DIR, 'analysis.html'), 'r') as f:
            template = Template(f.read())

        html_figure = f"""
        <div class="grid">
            <figure>
                <h3>Correlation of negative log E-value (nlogE) and normalized percentage of identity (NPID)</h3>
                <img src="{figures[-1][4]}" alt="Scatterplot and histogram of nlogE and NPID">
            </figure>
            {'<figure><h3>Maximum cluster size at different thresholds</h3><img src="' + figures[-1][5] + '" alt="Barplot of max cluster size"></figure>' if len(figures[-1]) == 6 else ''}
        </div>\n
        """
        self.report += html_figure
        
        for figs in figures:
            threshold = figs[0]
            fig_path_1 = figs[1]
            fig_path_2 = figs[2]

            fig_report = template.substitute(threshold=threshold,
                                            path_to_figure_1=fig_path_1,
                                            path_to_figure_2=fig_path_2)
            self.report += fig_report


    def save_html(self, output_file):
        """
        Parameters
        ----------
        output_file : str
            Output file
        """

        with open(os.path.join(HTML_DIR, 'page.html'), 'r') as f:
            template = Template(f.read())

        final_report = template.substitute(content=self.report)

        with open(output_file, 'w') as f:
            f.write(final_report)

    
        



