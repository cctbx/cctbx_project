import libtbx.load_env
from libtbx import test_utils

tst_list_base = []
try:
  import qttbx.qt  # noqa: F401
  tst_list_base.append("$D/regression/tst_phil_widgets.py")
  tst_list_base.append("$D/regression/tst_data_manager_widget.py")
except ImportError:
  pass

# Chat foundation tests — no Qt required. Alphabetized within the
# block so the list visibly drifts when a new test is added without
# updating this manifest.
tst_list_base.append("$D/regression/tst_chat_agent_core.py")
tst_list_base.append("$D/regression/tst_chat_diagnostics.py")
tst_list_base.append("$D/regression/tst_chat_logging_setup.py")
tst_list_base.append("$D/regression/tst_chat_mcp_client.py")
tst_list_base.append("$D/regression/tst_chat_profile.py")
tst_list_base.append("$D/regression/tst_chat_session.py")
tst_list_base.append("$D/regression/tst_chat_skills.py")
tst_list_base.append("$D/regression/tst_chat_storage.py")
tst_list_base.append("$D/regression/tst_chat_storage_session.py")
tst_list_base.append("$D/regression/tst_chat_tools.py")

# Qt-requiring chat widgets
try:
  import qttbx.qt  # noqa: F401
  tst_list_base.append("$D/regression/tst_chat_artifact_panel.py")
  tst_list_base.append("$D/regression/tst_chat_auto_height.py")
  tst_list_base.append("$D/regression/tst_chat_conversation_widgets.py")
  tst_list_base.append("$D/regression/tst_chat_credentials_dialog.py")
  tst_list_base.append("$D/regression/tst_chat_image.py")
  tst_list_base.append("$D/regression/tst_chat_markdown.py")
  tst_list_base.append("$D/regression/tst_chat_message_bubble.py")
  tst_list_base.append("$D/regression/tst_chat_message_input.py")
  tst_list_base.append("$D/regression/tst_chat_question_card.py")
  tst_list_base.append("$D/regression/tst_chat_runner.py")
  tst_list_base.append("$D/regression/tst_chat_tool_approval.py")
  tst_list_base.append("$D/regression/tst_chat_window_chrome.py")
except ImportError:
  pass

tst_list = tst_list_base

def run():
  build_dir = libtbx.env.under_build("qttbx")
  dist_dir = libtbx.env.dist_path("qttbx")
  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
