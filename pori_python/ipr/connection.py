import requests

import json
import os
import time
import zlib
from typing import Dict, List

from .constants import DEFAULT_URL
from .util import logger

IMAGE_MAX = 20  # cannot upload more than 20 images at a time


class IprConnection:
    def __init__(
        self,
        username: str,
        password: str,
        url: str = os.environ.get("IPR_URL", DEFAULT_URL),
    ):
        self.token = None
        self.url = url
        self.username = username
        self.password = password
        self.headers = {
            "Accept": "application/json",
            "Content-Type": "application/json",
            "Content-Encoding": "deflate",
        }
        self.cache: Dict[str, List[Dict]] = {}
        self.request_count = 0

    def request(self, endpoint: str, method: str = "GET", **kwargs) -> Dict:
        """Request wrapper to handle adding common headers and logging

        Args:
            endpoint (string): api endpoint, excluding the base uri
            method (str, optional): the http method. Defaults to 'GET'.

        Returns:
            dict: the json response as a python dict
        """
        url = f"{self.url}/{endpoint}"
        self.request_count += 1
        kwargs_header = kwargs.pop("headers", None)
        if kwargs_header:
            headers = json.loads(kwargs_header)
        else:
            headers = self.headers
        resp = requests.request(
            method, url, headers=headers, auth=(self.username, self.password), **kwargs
        )
        try:
            resp.raise_for_status()
        except requests.exceptions.HTTPError as err:
            # try to get more error details
            message = str(err)
            try:
                message += " " + resp.json()["error"]["message"]
            except Exception:
                pass

            raise requests.exceptions.HTTPError(message)
        if resp.status_code == 204:  # TODO: address this in api
            return {"status_code": 204}
        return resp.json()

    def post(self, uri: str, data: Dict = {}, **kwargs) -> Dict:
        """Convenience method for making post requests"""
        return self.request(
            uri,
            method="POST",
            data=zlib.compress(json.dumps(data, allow_nan=False).encode("utf-8")),
            **kwargs,
        )

    def get(self, uri: str, data: Dict = {}, **kwargs) -> Dict:
        """Convenience method for making get requests"""
        return self.request(
            uri,
            method="GET",
            data=zlib.compress(json.dumps(data, allow_nan=False).encode("utf-8")),
            **kwargs,
        )

    def delete(self, uri: str, data: Dict = {}, **kwargs) -> Dict:
        """Convenience method for making delete requests"""
        return self.request(
            uri,
            method="DELETE",
            data=zlib.compress(json.dumps(data, allow_nan=False).encode("utf-8")),
            headers=json.dumps({"Accept": "*/*"}),
            **kwargs,
        )

    def upload_report(
        self, content: Dict, mins_to_wait: int = 5, async_upload: bool = False
    ) -> Dict:
        if async_upload:
            # if async is used, the response for reports-async contains either 'jobStatus'
            # or 'report'. jobStatus is no longer available once the report is successfully
            # uploaded.
            initial_result = self.post("reports-async", content)

            report_id = initial_result["ident"]

            def check_status_result(result):
                if result.get("report", False):
                    return "upload complete"
                if result.get("jobStatus", False) and result["jobStatus"].get("state", False):
                    return result["jobStatus"]["state"]
                raise Exception(
                    "async report get returned with no report or jobStatus, or unexpected jobStatus type"
                )

            def check_status(interval: int = 5, num_attempts: int = 5):
                for i in range(num_attempts):
                    logger.info(f"checking report loading status in {interval} seconds")
                    time.sleep(interval)
                    current_status = self.get(f"reports-async/{report_id}")

                    check_result = check_status_result(current_status)

                    if check_result == "upload complete":
                        return current_status

                    if check_result == "failed":
                        raise Exception(
                            f"async report upload failed with reason: {current_status.get('jobStatus', {}).get('failedReason', 'Unknown')}"
                        )

                    if check_result not in [
                        "active",
                        "ready",
                        "waiting",
                        "completed",
                    ]:
                        raise Exception(f"async report upload in unexpected state: {check_result}")

                return current_status

            current_status = check_status()
            check_result = check_status_result(current_status)

            if check_result in ["active", "waiting"]:
                current_status = check_status(interval=30)
                check_result = check_status_result(current_status)

            if check_result in ["active", "waiting"]:
                current_status = check_status(interval=60, num_attempts=mins_to_wait)
                check_result = check_status_result(current_status)

            if check_result in ["active", "waiting"]:
                raise Exception(
                    f"async report upload taking longer than expected: {current_status}"
                )

            return current_status
        else:
            return self.post("reports", content)

    def set_analyst_comments(self, report_id: str, data: Dict) -> Dict:
        """
        Update report comments to an existing report

        TODO:
            Add to main upload.
            Pending: https://www.bcgsc.ca/jira/browse/DEVSU-1177
        """
        return self.request(
            f"/reports/{report_id}/summary/analyst-comments",
            method="PUT",
            data=zlib.compress(json.dumps(data, allow_nan=False).encode("utf-8")),
        )

    def post_images(self, report_id: str, files: Dict[str, str], data: Dict[str, str] = {}) -> None:
        """
        Post images to the report
        """
        file_keys = list(files.keys())
        start_index = 0
        image_errors = set()
        while start_index < len(file_keys):
            current_files = {}
            for key in file_keys[start_index : start_index + IMAGE_MAX]:
                path = files[key]
                if not os.path.exists(path):
                    raise FileNotFoundError(path)
                current_files[key] = path
            open_files = {k: open(f, "rb") for (k, f) in current_files.items()}
            try:
                resp = self.request(
                    f"reports/{report_id}/image",
                    method="POST",
                    data=data,
                    files=open_files,
                    headers=json.dumps({}),
                )
                for status in resp:
                    if status.get("upload") != "successful":
                        image_errors.add(status["key"])
            finally:
                for handler in open_files.values():
                    handler.close()
            start_index += IMAGE_MAX
        if image_errors:
            raise ValueError(f'Error uploading images ({", ".join(sorted(list(image_errors)))})')

    def get_spec(self) -> Dict:
        """
        Get the current IPR spec, for the purposes of current report upload fields
        """
        return self.request("/spec.json", method="GET")
