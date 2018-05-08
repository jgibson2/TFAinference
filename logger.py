import datetime


###
# Class Logger:
#	Logs to a file
###

class Logger:
    fh = None

    ###
    # Constructor:
    #	Args:
    #		filehandle : filehandle of file to log to
    #		truncate : Whether to truncate the file or not; defaults to False
    #	Returns:
    #		Object reference
    ###
    def __init__(self, filehandle, truncate=False):
        if type(filehandle) == str:
            try:
                filehandle = open(filehandle, 'w')
            except:
                raise ValueError('Unable to open file for writing!')
        self.fh = filehandle
        time = datetime.datetime.now()
        if truncate:
            self.fh.truncate()
        self.log("Log started")

    ###
    # log:
    #	Args:
    #		x : object to log.
    #	Returns:
    #		None
    #	Logs x to class filehandle
    ###
    def log(self, x):
        try:
            time = datetime.datetime.now()
            self.fh.write(time.strftime('%X: ' + str(x) + '\n'))
        except ValueError as err:
            time = datetime.datetime.now()
            try:
                self.fh.write(time.strftime('%X: ' + repr(x) + '\n'))
            except ValueError as internal_err:
                self.fh.write(time.strftime('%X: Could not write object to log.\n'))
