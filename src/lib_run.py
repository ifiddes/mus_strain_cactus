import os
import subprocess


def RunCommands(cmds, local_temp_dir, in_pipes=None, out_pipes=None,
                 err_pipes=None, ignore_returns=False, serial=False):
  """ Wrapper for RunCommandsS and RunCommandsP.
  """
  if serial:
    RunCommandsS(cmds, local_temp_dir, in_ipes, out_pipes, err_pipes,
                 ignore_returns)
  else:
    RunCommandsP(cmds, local_temp_dir, in_ipes, out_pipes, err_pipes,
                 ignore_returns)


def HandlePipes(in_pipes, out_pipes, err_pipes, ignore_returns):
  """ Generate the correct data structures for RunCommands*() functions.
  """
  if in_pipes is None:
    in_pipes = [None] * len(cmds)
  if out_pipes is None:
    out_pipes = [None] * len(cmds)
  if err_pipes is None:
    err_pipes = [None] * len(cmds)
  if ignore_returns is None:
    ignore_returns = [False] * len(cmds)
  return in_pipes, out_pipes, err_pipes, ignore_returns


def HandlePipesInstance(in_pipe, out_pipe, err_pipe):
  """ Generate the correct data structures for a single command.
  """
  if in_pipe is None:
    sin = None
  else:
    sin = subprocess.PIPE
  if out_pipe is None:
    sout = None
  else:
    sout = subprocess.PIPE
  if err_pipe is None:
    serr = None
  else:
    serr = subprocess.PIPE
  return sin, sout, serr


def RunCommandsS(cmds, local_temp_dir, in_pipes=None, out_pipes=None,
                 err_pipes=None, ignore_returns=None):
  """ Uses the subprocess module to issue serial processes from the cmds list.
  Arguments:
    ignore_returns: if true, return code is ignored
  """
  if not os.path.exists(local_temp_dir):
    raise ValueError('local_temp_dir "%s" does not exist.' % local_temp_dir)
  in_pipes, out_pipes, err_pipes, ignore_returns = HandlePipes(
    in_pipes, out_pipes, err_pipes, ignore_returns)
  for i, c in enumerate(cmds, 0):
    sin, sout, serr = HandlePipesInstance(in_pipes[i], out_pipes[i],
                                          err_pipes[i])
    p = subprocess.Popen(c, cwd=local_temp_dir, stdin=sin,
                         stdout=sout, stderr=serr)
    if in_pipes[i] is None:
      sin = None
    else:
      if not os.path.exists(in_pipes[i]):
        raise IOError('Unable to locate inPipe file: %s for command %s'
                      % (in_pipes[i], ' '.join(c)))
      sin = open(in_pipes[i], 'r').read()
    if out_pipes[i] is None:
      pout, perr = p.communicate(sin)
      if not ignore_returns[i]:
        HandleReturnCode(p.returncode, cmds[i])
    else:
      f = open(out_pipes[i], 'w')
      f.write(p.communicate(sin)[0])
      f.close()
      if not ignore_returns[i]:
        HandleReturnCode(p.returncode, cmds[i])


def RunCommandsP(cmds, local_temp_dir, in_pipes=None, out_pipes=None,
                 err_pipes=None, ignore_returns=None, **kwargs):
  """ Uses the subprocess module to issue parallel processes from cmds list.
  """
  if not os.path.exists(local_temp_dir):
    raise ValueError('local_temp_dir "%s" does not exist.' % local_temp_dir)
  procs = []
  in_pipes, out_pipes, err_pipes, ignore_returns = HandlePipes(
    in_pipes, out_pipes, err_pipes, ingore_returns)
  for i, c in enumerate(cmds, 0):
    sin, sout, serr = HandlePipesInstance(in_pipes[i], out_pipes[i],
                                          err_pipes[i])
    procs.append(subprocess.Popen(
        c, cwd=local_temp_dir, stdin=sin, stdout=sout, stderr=serr))
  for i, p in enumerate(procs, 0):
    if in_pipes[i] is None:
      sin = None
    else:
      if not os.path.exists(in_pipes[i]):
        raise IOError('Unable to locate inPipe file: %s for command %s'
                      % (in_pipes[i], cmds[i]))
      sin = open(in_pipes[i], 'r').read()
    if out_pipes[i] is None:
      pout, perr = p.communicate(sin)
      if not ignore_returns[i]:
        HandleReturnCode(p.returncode, cmds[i])
    else:
      f = open(out_pipes[i], 'w')
      f.write(p.communicate(sin)[0])
      f.close()
      if not ignore_returns[i]:
        HandleReturnCode(p.returncode, cmds[i] + ['< %s > %s 2> %s'
                                                  % (in_pipes[i],
                                                     out_pipes[i],
                                                     err_pipes[i])])


def HandleReturnCode(retcode, cmd):
  """ Handle the return codes from RunCommands.
  """
  if not isinstance(retcode, int):
    raise TypeError('HandleReturnCode takes an integer for '
                    'retcode, not a %s.' % retcode.__class__)
  if not isinstance(cmd, list):
    raise TypeError('HandleReturnCode takes a list for '
                    'cmd, not a %s.' % cmd.__class__)
  if retcode:
    if retcode < 0:
      raise RuntimeError('Experienced an error while trying to execute: '
                         '%s SIGNAL:%d' %(' '.join(cmd), -retcode))
    else:
      raise RuntimeError('Experienced an error while trying to execute: '
                         '%s retcode:%d' %(' '.join(cmd), retcode))


def Which(program, extra_path_list=None):
  """ Which() acts like the unix utility which, but is portable between os.
  If the program does not exist in the PATH then 'None' is returned.
  extra_path_list is an os-appropritate list of paths to start the search
  for the executable.
  """
  if extra_path_list is not None:
    if not isinstance(extra_path_list, list):
      raise RuntimeError('extra_path_list must be a list, not a %s'
                         % extra_path_list.__class__)
  def is_exe(fpath):
    return os.path.exists(fpath) and os.access(fpath, os.X_OK)
  fpath, fname = os.path.split(program)
  if fpath != '':
    if is_exe(program):
      return program
  else:
    if extra_path_list is not None:
      for path in extra_path_list:
        exe_file = os.path.join(path, program)
        if is_exe(exe_file):
          return exe_file
    for path in os.environ["PATH"].split(os.pathsep):
      exe_file = os.path.join(path, program)
      if is_exe(exe_file):
        return exe_file
    exe_file = os.path.join(os.getcwd(), program)
    if is_exe(exe_file):
      return exe_file
  return None
