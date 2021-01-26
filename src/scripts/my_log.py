import sys
import os
import logging

# if not os.path.exists('log'):
#     os.makedirs('log')

# By default the root logger is set to WARNING and all loggers you define inherit that value.
logging.root.setLevel(logging.NOTSET)

logger = logging.getLogger(__name__)

c_handler = logging.StreamHandler()
f_handler = logging.FileHandler('run.log')

f_format = logging.Formatter('%(asctime)s - [%(levelname)s] %(filename)s > (Method: %(funcName)s, line no. %(lineno)d) :: %(message)s')
f_handler.setFormatter(f_format)

c_handler.setLevel(logging.DEBUG)
f_handler.setLevel(logging.DEBUG)

logger.addHandler(c_handler)
logger.addHandler(f_handler)