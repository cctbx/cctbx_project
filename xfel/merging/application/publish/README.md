# Publish mtz and log files

### Synopsis

This worker publishes merging results (mtz and main log) to a Google Drive folder. Subfolders are created for
the dataset name and version (e.g. `dataset1/v001`). Authentication is handled via Google Cloud service accounts,
which can be created easily with limited permissions, so that the credential file for the uploading account can be
shared safely.

The rest of this file details the setup of a service account and a sharing destination folder. The steps may depend
on the Google services your organization subscribes to; they have been tested for an LBNL account.

### Setup

1. Create a Drive folder as a destination for MTZ files. **NOTE 2025/7/11**: This procedure has changed to reflect
   new rules for storage quotas.
   1. The upload destination must be a Shared Drive, which is different from a regular drive that has been shared
      with other accounts. Since most experiments use the second (incorrect) model, we will create a separate destination
      in a Shared Drive and link it inside the experiment folder.
   2. Visit https://drive.google.com/drive/shared-drives and create a new Shared Drive, which may be reused for
      multiple experiments (example name: `Uploads`). Inside this Shared Drive, create a folder which will be the
      destination for the current experiment (example name: `mfx101080424_mtz`).
   3. Right-click the new folder and select Organize-->Add shortcut. Place the shortcut inside the experiment folder.
      Make sure to update the sharing permissions on this folder to match the experiment folder.
   4. Note the identifier string of the destination folder. If its URL is
      https://drive.google.com/drive/u/0/folders/16-VdcQw_9Jy9MxDDdjTxpYX5iiE7nsk, then its ID string is 
      16-VdcQw_9Jy9MxDDdjTxpYX5iiE7nsk.
2. Create a Google Cloud project. Try visiting: https://console.cloud.google.com/cloud-resource-manager
and clicking “Create Project”. A convenient name would be e.g. `dwpaley-mtz-upload`. You will be able to reuse this
project to create service accounts for multiple beamtimes.
3. Enable the Google Drive API for the project.
   1. Visit: https://console.cloud.google.com/apis/dashboard. Select the project. Click “Enable APIs and Services”,
then search for Drive. Click the Enable button.
4. Create a service account in the project.
   1. Try visiting: https://console.cloud.google.com/iam-admin/serviceaccounts. Select the project you just created.
Click “CREATE SERVICE ACCOUNT”. Give the account a name connected to the beamtime, e.g. `mfxlv4318-uploader`.
   2. Share the upload destination folder (from step 1) with the service account: Note the email address for the
service account. Visit the `mfx101080424_mtz` Drive folder and share it with this email address.
   3. While you are still on the page “IAM & Admin-->Service Accounts”, click the “Actions” button for the new service
account and select “Manage Keys”. Generate a new key by clicking “Add Key”, “Create new key”. Save the resulting
JSON file and transfer it to a convenient location on the cluster where the data will be processed.
5. Add `publish` as the final step in the merging phil item `dispatch.step_list`. Add the merging phil parameters
`publish.drive.credential_file` (from step 4.iii) and `publish.drive.shared_folder_id` (from step 1).
