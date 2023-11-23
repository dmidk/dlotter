import os
from dmit import dmitio

class s3util:
    """Class for uploading files to s3"""

    def __init__(self) -> None:
        """Constructor for s3 class
        """

        self.AWS_S3_BUCKET = os.environ.get('AWS_S3_BUCKET')

        return


    def send_files_to_s3(self, files:list, s3path:str) -> None:
        """Send files to s3

        Parameters
        ----------
        files : list
            List of files to send to s3
        s3path : str
            Path to s3 bucket

        Returns
        -------
        None
            None
        """

        S3 = dmitio.s3(cert=False)

        # Sort the files to make sure they are uploaded in the correct order
        files = list(sorted(files))

        # TODO: # If more than one analysis is present, this wont work
        k=0
        for file in files:

            # Use split to extract the date from the filename
            date = file.split('_')[1]
            time = file.split('_')[2].split('-')[0]
            analysis = date + time

            filename = os.path.basename(file).split('_')[0]+'_'+str(k)+'.png'

            remote_file = os.path.join(s3path, analysis, filename)

            S3.upload(self.AWS_S3_BUCKET, remote_file, file)
            k+=1

        return