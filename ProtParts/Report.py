
from string import Template
import os

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

        with open('template/params.html', 'r') as f:
            template = Template(f.read())
        
        params_report = template.substitute(threshold_c=params.threshold_c,
                                            threshold_r=params.threshold_r,
                                            num_partitions=params.num_partitions,
                                            redundancy_reduction=bool(params.threshold_r),
                                            cluster_separation=params.separation)
        
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

        
        # Add a row to the HTML table for each row in the CSV file
        for i, row in enumerate(results):
            if i == 0:  # Assuming the first row is the header
                html_table += '  <thead>\n    <tr>' + ''.join([f'<th>{cell}</th>' for cell in row]) + '</tr>\n  </thead>\n  <tbody>\n'
            else:
                html_table += '    <tr>' + ''.join([f'<td>{cell}</td>' for cell in row[:-1]]) + f'<td><a href="{row[-1]}" download>Link</a></td></tr>\n'
        html_table += '  </tbody>\n</table>\n'


        with open('template/results.html', 'r') as f:
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

        self.report += "<h2>Analysis</h2>\n"

        # load the template
        with open('template/analysis.html', 'r') as f:
            template = Template(f.read())

        if len(figures[0]) == 5:
            html_figure = f"""
            <div class="grid">
            <figure>
            <h3>Maximum cluster size at different thresholds</h3>
            <img src="{figures[0][4]}" alt="Barplot of max cluster size">
            </figure>
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

        with open('template/page.html', 'r') as f:
            template = Template(f.read())

        final_report = template.substitute(content=self.report)

        with open(output_file, 'w') as f:
            f.write(final_report)

    
        



