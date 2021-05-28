#!/usr/bin/env python
"""Runs W-Phase, sends results to seiscomp3 messaging and copies files to S3."""
from __future__ import absolute_import, print_function

import logging
import os
import sys
import json

from datetime import datetime, timedelta

from obspy.core import UTCDateTime

from seiscomp3 import DataModel as DM, Logging, IO, Core
from seiscomp3.Client import Application

from wphase import logger, runwphase, settings
from wphase.email import send_email
from wphase.result import WPhaseParser, FMItem, jsonencode_np
from wphase.seiscomp import createAndSendObjects, writeSCML


# this is used in determining the name of the log file that this application
# logs to. I think this must be the name of this file, but cannot be sure as the
# logging is Gempa's black magic.
APPLICATION_NAME = 'wphase'

LOG_FILE_NAME = '{}/.seiscomp3/log/{}.log'.format(
    os.environ['HOME'], APPLICATION_NAME)

class LogRelay(logging.Handler):
    """Handler that forwards a python logger to seiscomp logs"""
    def emit(self, record):
        msg = self.format(record)
        if record.levelname == 'DEBUG':
            Logging.debug(msg)
        elif record.levelname == 'INFO':
            Logging.info(msg)
        elif record.levelname == 'WARNING':
            Logging.warning(msg)
        else:
            Logging.error(msg)

# Send W-Phase logs and python warnings to seiscomp logs
logger.setLevel(logging.DEBUG)
logger.addHandler(LogRelay(level=logging.DEBUG))
wlogger = logging.getLogger('py.warnings')
wlogger.setLevel(logging.WARNING)
wlogger.addHandler(LogRelay(level=logging.WARNING))
logging.captureWarnings(True)

class WPhase(Application):
    def __init__(self, argc, argv):
        # Default to loglevel 1
        Application.__init__(self, argc, argv)

        # default location to write outputs to
        self.output = settings.DEFAULT_OUTPUT_DIR
        self.filename = None
        self.mag_type = None
        self.mag_value = None
        self.server = 'IRIS'
        self.networks = 'ALL'
        self.region = 'not specified'
        self.evid = None
        self.resultid = None
        self.notificationemail = None
        self.fromemail = None
        self.email_aws_region = None
        self.email_method = 'ses'

        self.email_subject_postfix = ''
        self.email_subject_prefix = ''

        self.smtp_server = None
        self.smtp_port = 25
        self.smtp_ssl = False
        self.smtp_tls = False
        self.smtp_user = None
        self.smtp_password = None

        self.write_s3 = False
        self.bucket_name = None
        self.agency = 'GA'
        self.make_maps = True
        self.overwrite = False
        self.save_waveforms = None
        self.save_inventory = None
        self.waveforms = None
        self.inventory = None

        # enable messaging support
        self.setMessagingEnabled(True)

        # disable database access
        self.setDatabaseEnabled(False, False)

        # default spread username
        self.setMessagingUsername("gawphase")

        # send all objects to the focal mechanism group
        self.setPrimaryMessagingGroup("FOCMECH")



    def createCommandLineDescription(self):
        """
        Override this method (including a super call) in a subclass to add
        custom options.
        """

        self.commandline().addGroup("Output")
        self.commandline().addStringOption(
            "Output",
            "outputs,o",
            "Directory to write output to. Defaults to /tmp/wphase-output.")
        self.commandline().addStringOption(
            "Output",
            "save-waveforms",
            "Path in which to save raw waveforms (as miniseed)")
        self.commandline().addStringOption(
            "Output",
            "save-inventory",
            "Path in which to save raw inventory (as stationXML)")

        self.commandline().addGroup("Input")
        self.commandline().addStringOption(
            "Input",
            "lat",
            "Latitude of the event.")
        self.commandline().addStringOption(
            "Input",
            "lon",
            "Longitude of the event.")
        self.commandline().addStringOption(
            "Input",
            "depth",
            "Depth of the event.")
        self.commandline().addStringOption(
            "Input",
            "sourcezone",
            "Source zone of the event.")
        self.commandline().addStringOption(
            "Input",
            "time",
            "Time of the event.")
        self.commandline().addStringOption(
            "Input",
            "magtype",
            "The type of magnitude of the triggering event.")
        self.commandline().addStringOption(
            "Input",
            "magvalue",
            "The magnitude of the triggering event.")
        self.commandline().addStringOption(
            "Input",
            "waveforms",
            "Path to waveforms to use (instead of FDSN)")
        self.commandline().addStringOption(
            "Input",
            "inventory",
            "Path to inventory to use (instead of FDSN)")
        self.commandline().addStringOption(
            "Input",
            "server",
            "The FDSN server to use. This can be any server supported by " + \
            "the obspy FDSN client, or an arbitrar URL.")
        self.commandline().addStringOption(
            "Input",
            "networks",
            "A comma separated list of networks to use.")
        self.commandline().addStringOption(
            "Input",
            "region",
            "The name of the region in which the event occured.")
        self.commandline().addStringOption(
            "Input",
            "evid,E",
            "Event ID.")
        self.commandline().addStringOption(
            "Input",
            "resultid",
            "A string that (along with the evid) will uniquely identify this particular w-phase result. ")
        self.commandline().addStringOption(
            "Input",
            "notificationemail",
            "Email address to send notification to.")
        self.commandline().addStringOption(
            "Input",
            "fromemail",
            "Email address to send notification from.")
        self.commandline().addStringOption(
            "Input",
            "emailmethod",
            "Method to use to send notification email (either ses or smtp)")
        self.commandline().addStringOption(
            "Input",
            "emailawsregion",
            "If using SES, AWS Region to send email notifications from.")
        self.commandline().addStringOption(
            "Input",
            "email-subject-prefix",
            "String to add to start of notification email subject.")
        self.commandline().addStringOption(
            "Input",
            "email-subject-postfix",
            "String to add to end of notification email subject.")
        self.commandline().addStringOption(
            "Input",
            "smtp-server",
            "SMTP server to use to send notification email.")
        self.commandline().addStringOption(
            "Input",
            "smtp-port",
            "SMTP port.")
        self.commandline().addStringOption(
            "Input",
            "smtp-user",
            "SMTP username.")
        self.commandline().addStringOption(
            "Input",
            "smtp-password",
            "SMTP username.")
        self.commandline().addOption(
            "Input",
            "smtp-ssl",
            "Enable SSL for SMTP.")
        self.commandline().addOption(
            "Input",
            "smtp-tls",
            "Enable TLS for SMTP.")
        self.commandline().addStringOption(
            "Input",
            "writeS3",
            "Should results be written to S3?")
        self.commandline().addStringOption(
            "Input",
            "bucketname",
            "The name of the S3 bucket to write to.")
        self.commandline().addStringOption(
            "Input",
            "agency",
            "Agency code for SC3 output")
        self.commandline().addOption(
            "Input",
            "make-maps",
            "Whether to create map plots or not.")
        self.commandline().addOption(
            "Input",
            "overwrite",
            "Whether to overwrite existing outputs or not.")

    def processArg(self, name, to=None, default=None, conv=str):
        if to is None:
            to = name

        try:
            val = self.commandline().optionString(name)
        except Exception:
            if default is not None:
                setattr(self, to, default)
        else:
            setattr(self, to, conv(val))


    def validateParameters(self):
        """Called by the seiscomp client Application setup."""
        if not Application.validateParameters(self):
            return False

        getter = self.processArg
        def getflag(name, to=None):
            getter(name, to=to, conv=lambda x: True)


        try:
            # If there is an unrecognized option it must be a JSON file
            # of wphase outputs. In this case, the file is parsed and pushed
            # to the messaging system and written to disk.
            self.filename = self.commandline().unrecognizedOptions()[0]
        except Exception:
            # Otherwise we expect a description of the location. Wphase is
            # then run and results pushed pushed to the messaging system
            # and written to disk.

            # depth only needed for BoM XML
            try: depth = float(self.commandline().optionString("depth"))
            except Exception: depth = 0.

            getter('save-waveforms', 'save_waveforms')
            getter('save-inventory', 'save_inventory')
            getter('waveforms', 'waveforms')
            getter('inventory', 'inventory')

            try:
                self.eqinfo = {
                    'lat' : float(self.commandline().optionString("lat")),
                    'lon' : float(self.commandline().optionString("lon")),
                    'dep' : depth,
                    'time': UTCDateTime(self.commandline().optionString("time"))}
            except Exception:
                logger.error('You must provide either lat/lon/time or a JSON payload')
                return False

            getter('sourcezone')
            getter('magtype', 'mag_type')
            getter('magvalue', 'mag_value', conv=float)
            getter('outputs', 'output')
            getter('server')
            getter('networks')
            getter('region')
            getter('evid')
            getter('resultid')
            getter('notificationemail')
            getter('emailmethod', 'email_method')
            getter('fromemail')
            getter('email-subject-prefix', 'email_subject_prefix')
            getter('email-subject-postfix', 'email_subject_postfix')
            getter('emailawsregion', 'email_aws_region', 'us-west-2')
            getter('writeS3', 'write_s3', False, lambda x: True)
            getter('bucketname', 'bucket_name')
            getter('agency')
            getflag('make-maps', 'make_maps')
            getflag('overwrite', 'overwrite')

            getter('smtp-server', 'smtp_server')
            getter('smtp-port', 'smtp_port')
            getter('smtp-user', 'smtp_user')
            getter('smtp-password', 'smtp_password')
            getflag('smtp-ssl', 'smtp_ssl')
            getflag('smtp-tls', 'smtp_tls')

            if self.evid is not None:
                self.output = os.path.join(self.output, self.evid)

            if self.resultid is not None:
                self.output = os.path.join(self.output, self.resultid)

            if self.notificationemail is not None:
                self.notificationemail = str(self.notificationemail).split(',')

            if self.write_s3 and (
                    self.evid is None or \
                    self.bucket_name is None):
                logger.error('attempt to write to s3, but did not provide bucket name.')
                return False

            if self.notificationemail is not None and (
                    self.mag_type is None or \
                    self.mag_value is None or \
                    self.evid is None or \
                    self.resultid is None
                    ):
                logger.error('cannot send email.')
                return False

        return True



    def init(self):
        """
        SC3 specific method.

        Returning False means that we do not enter the SeisComP3 run loop.
        """

        if Application.init(self) == False:
            return False

        item = None
        res = None
        parser = WPhaseParser()

        Logging.enableConsoleLogging(Logging.getGlobalChannel("error"))
        wphase_failed = False

        if self.evid is not None:
            self.eqinfo['id'] = self.evid
        if self.filename is None:
            try:
                logger.info("Starting W-Phase.")
                res = runwphase(
                    output_dir=self.output,
                    server=self.server,
                    eqinfo=self.eqinfo,
                    networks=self.networks,
                    make_maps=self.make_maps,
                    output_dir_can_exist=self.overwrite,
                    save_waveforms=self.save_waveforms,
                    save_inventory=self.save_inventory,
                    waveforms=self.waveforms,
                    inventory=self.inventory)
            except Exception:
                from traceback import format_exc
                logger.error('failed to run wphase: %s', format_exc())
                wphase_failed = True
            else:
                if self.evid is not None:
                    try:
                        # TODO: Should this be done in runwphase?
                        res[settings.WPHASE_EVENT_KEY]['id'] = self.evid
                        with open(os.path.join(
                                self.output,
                                settings.WPHASE_OUTPUT_FILE_NAME), 'w') as of:
                            json.dump(res, of, default=jsonencode_np)
                    except Exception as e:
                        # not sure how we would get here, but we just don't want
                        # to stop the rest of processing
                        logger.error('failed to add event id to event: %s', e)

                try:
                    try: res_dict = res.as_dict()
                    except Exception: res_dict = res
                    item = parser.read(json_data=res_dict)
                except Exception as e:
                    logger.error('failed to parse event JSON for SC3: %s', e)
                    item = None

        else:
            try:
                item = parser.read(filename=self.filename)
            except Exception as e:
                logger.error('failed parse event JSON for SC3: %s', e)

        if item is not None:
            try:
                objs = createAndSendObjects(item, self.connection(),
                                            evid=self.evid, logging=logger,
                                            agency=self.agency)
            except Exception as e:
                logger.error('failed create objects for SC3: %s', e)
            else:
                try:
                    # write output to file
                    filename = os.path.join(self.output, 'sc3.xml')
                    writeSCML(filename, objs)
                    logger.info("Stored results in SCML as %s", filename)
                except Exception as e:
                    logger.error('failed write SCML to file: %s', e)


        if self.write_s3:
            # We always want to write to S3, even in the event of failure as
            # the JSON output may explain the failure. We also want to do it
            # at the very end since we keep the sc3 log file.

            try:
                from wphase.aws import write_to_s3
                write_to_s3(
                    self.output,
                    self.bucket_name,
                    self.evid,
                    self.resultid,
                    [(LOG_FILE_NAME, 'sc3.log')],
                    logger.error)
            except Exception as e:
                logger.error('failed write to S3: %s', e)
            finally:
                # since we may do more runs, remove the log so we don't
                # contaminate later runs with info for this run
                try:
                    os.remove(LOG_FILE_NAME)
                except OSError:
                    # ... don't log this
                    pass

        if self.notificationemail:
            # must be done after writing to S3
            success = res is not None and 'MomentTensor' in res
            subject, body = self.createEmail(event_id=self.evid,
                                             result_id=self.resultid,
                                             result_dict=res,
                                             call_succeeded=success,
                                             )
            send_email(recipients=self.notificationemail,
                       subject=self.email_subject_prefix + str(subject)
                               + self.email_subject_postfix,
                       message=body,
                       from_email=self.fromemail,
                       method=self.email_method,

                       email_aws_region=self.email_aws_region,

                       server=self.smtp_server,
                       port=self.smtp_port,
                       user=self.smtp_user,
                       password=self.smtp_password,
                       ssl=self.smtp_ssl,
                       tls=self.smtp_tls,
                       )

        sys.exit(1 if wphase_failed else 0)

    def createEmail(self, event_id, result_id, result_dict, call_succeeded):
        """Create email subject and body to notify result.

        Override this in a subclass to produce a custom email."""
        subject = 'W-Phase succeeded' if call_succeeded else 'W-Phase failed'
        body = "<h2>W-Phase inversion {} for {} results:</h2> <pre>{}</pre>".format(result_id, event_id, json.dumps(result_dict))
        return subject, body


def main():
    argv = sys.argv
    app = WPhase(len(argv), argv)
    sys.exit(app())


if __name__ == "__main__":
    main()
