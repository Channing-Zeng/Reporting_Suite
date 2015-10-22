from os.path import join, dirname
from traceback import format_exc

from source.logger import info, err, warn
from source.tools_from_cnf import get_system_path


def connect_to_server(server_url, username, password):
    # html_report_url = 'http://ngs.usbod.astrazeneca.net/reports/' + bcbio_structure.project_name + '/' + \
    #     relpath(html_report_fpath, bcbio_structure.final_dirpath)
    info('Connecting to the server ' + server_url + '...')

    rsa_key_path = get_system_path(None, join(dirname(__file__), 'id_rsa'), is_critical=False)
    if not rsa_key_path:
        err('Could not find key ' + rsa_key_path)

    try:
        from ext_modules.paramiko import SSHClient, RSAKey, AutoAddPolicy
    except ImportError as e:
        err(format_exc())
        warn()
        err('Cannot improt SSHClient - skipping trasnferring symlinking to the ngs-website')
        warn()
    else:
        ssh = SSHClient()
        ssh.load_system_host_keys()
        # ki = RSAKey.from_private_key_file(filename=rsa_key_path)
        ssh.set_missing_host_key_policy(AutoAddPolicy())
        try:
            key = RSAKey(filename=rsa_key_path, password='%1!6vLaD')
        except:
            warn('Cannot read RSAKey from ' + rsa_key_path)
            warn()
            err(format_exc())
            warn()
        else:
            info('Succesfully read RSAKey from ' + rsa_key_path)
            try:
                ssh.connect(server_url, username=username, password=password, pkey=key)
            except:
                warn('Cannot connect to ' + server_url + ':')
                warn()
                err(format_exc())
                warn()
            else:
                info('Succesfully connected to ' + server_url)
                return ssh

    return None