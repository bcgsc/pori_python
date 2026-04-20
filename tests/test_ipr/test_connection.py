import os
import pytest
from unittest import mock

from pori_python.ipr.connection import IprConnection

IMAGE_DIR = os.path.join(os.path.dirname(__file__), '../../docs/images')


class TestPostImages:
    def test_no_images_ok(self):
        def request(*args, **kwargs):
            m = mock.MagicMock(
                json=lambda: [{'upload': 'successful'}], raise_for_status=lambda: None
            )
            return m

        with mock.patch('pori_python.ipr.connection.requests.request', request):
            conn = IprConnection('user', 'pass')
            result = conn.post_images('report_id', files={}, data={})
            assert result is None

    def test_images_load_ok(self):
        def request(*args, **kwargs):
            m = mock.MagicMock(
                json=lambda: [{'upload': 'successful'}], raise_for_status=lambda: None
            )
            return m

        with mock.patch('pori_python.ipr.connection.requests.request', request):
            conn = IprConnection('user', 'pass')
            result = conn.post_images(
                'report_id',
                files={
                    'expression.correlation': os.path.join(IMAGE_DIR, 'expression_correlation.png'),
                    'mixcr.circos_trb_vj_gene_usage': os.path.join(
                        IMAGE_DIR, 'mixcr.circos_trb_vj_gene_usage.png'
                    ),
                },
                data={},
            )
            assert result is None

    def test_images_with_data_load_ok(self):
        def request(*args, **kwargs):
            m = mock.MagicMock(
                json=lambda: [{'upload': 'successful'}], raise_for_status=lambda: None
            )
            return m

        with mock.patch('pori_python.ipr.connection.requests.request', request):
            conn = IprConnection('user', 'pass')
            result = conn.post_images(
                'report_id',
                files={
                    'expression.correlation': os.path.join(IMAGE_DIR, 'expression_correlation.png'),
                    'mixcr.circos_trb_vj_gene_usage': os.path.join(
                        IMAGE_DIR, 'mixcr.circos_trb_vj_gene_usage.png'
                    ),
                },
                data={'expression.correlation.title': 'this is a title'},
            )
            assert result is None

    def test_bad_file(self):
        def request(*args, **kwargs):
            m = mock.MagicMock(
                json=lambda: [{'upload': 'successful'}], raise_for_status=lambda: None
            )
            return m

        with mock.patch('pori_python.ipr.connection.requests.request', request):
            conn = IprConnection('user', 'pass')
            with pytest.raises(FileNotFoundError):
                conn.post_images(
                    'report_id', files={'expression.correlation': 'thing/that/does/not/exist.png'}
                )

    def test_failed_image_load(self):
        def request(*args, **kwargs):
            m = mock.MagicMock(
                json=lambda: [{'upload': 'anything else', 'key': 'thing'}],
                raise_for_status=lambda: None,
            )
            return m

        with mock.patch('pori_python.ipr.connection.requests.request', request):
            conn = IprConnection('user', 'pass')
            with pytest.raises(ValueError):
                conn.post_images(
                    'report_id',
                    {
                        'expression.correlation': os.path.join(
                            IMAGE_DIR, 'expression_correlation.png'
                        )
                    },
                )


class TestCheckUploadPermission:
    def _user_response(self, groups=None, projects=None):
        return {
            'groups': [{'name': g} for g in (groups or [])],
            'projects': [{'name': p} for p in (projects or [])],
        }

    def test_rejects_user_without_create_report_access(self):
        conn = IprConnection('user', 'pass')
        conn.get = mock.MagicMock(
            side_effect=[[{'name': 'TEST'}], self._user_response(projects=['TEST'])]
        )
        conn.post = mock.MagicMock()

        with pytest.raises(Exception, match='User does not have report creation permission'):
            conn.check_upload_permission('TEST')

        conn.post.assert_not_called()

    def test_rejects_user_without_project_access(self):
        conn = IprConnection('user', 'pass')
        conn.get = mock.MagicMock(
            side_effect=[
                [{'name': 'TEST'}],
                self._user_response(groups=['create report access'], projects=['OTHER']),
            ]
        )
        conn.post = mock.MagicMock()

        with pytest.raises(Exception, match='User has no permission to create report in project TEST'):
            conn.check_upload_permission('TEST')

        conn.post.assert_not_called()

    def test_allows_user_with_project_and_create_report_access(self):
        conn = IprConnection('user', 'pass')
        conn.get = mock.MagicMock(
            side_effect=[
                [{'name': 'TEST'}],
                self._user_response(groups=['create report access'], projects=['TEST']),
            ]
        )
        conn.post = mock.MagicMock()

        conn.check_upload_permission('TEST')

        conn.post.assert_not_called()

    def test_manager_has_implicit_create_report_access(self):
        conn = IprConnection('user', 'pass')
        conn.get = mock.MagicMock(
            side_effect=[
                [{'name': 'TEST'}],
                self._user_response(groups=['manager'], projects=['TEST']),
            ]
        )
        conn.post = mock.MagicMock()

        conn.check_upload_permission('TEST')

        conn.post.assert_not_called()

    def test_admin_bypasses_all_checks(self):
        conn = IprConnection('user', 'pass')
        conn.get = mock.MagicMock(
            side_effect=[
                [{'name': 'TEST'}],
                self._user_response(groups=['admin'], projects=[]),
            ]
        )
        conn.post = mock.MagicMock()

        conn.check_upload_permission('TEST')

        conn.post.assert_not_called()

    def test_admin_creates_missing_project(self):
        conn = IprConnection('user', 'pass')
        conn.get = mock.MagicMock(
            side_effect=[
                [{'name': 'OTHER'}],
                self._user_response(groups=['admin'], projects=[]),
            ]
        )
        conn.post = mock.MagicMock()

        conn.check_upload_permission('TEST')

        conn.post.assert_called_once_with('project', {'name': 'TEST'})

    def test_all_projects_access_without_project_membership(self):
        conn = IprConnection('user', 'pass')
        conn.get = mock.MagicMock(
            side_effect=[
                [{'name': 'TEST'}],
                self._user_response(
                    groups=['create report access', 'all projects access'], projects=[]
                ),
            ]
        )
        conn.post = mock.MagicMock()

        conn.check_upload_permission('TEST')

        conn.post.assert_not_called()

    def test_creates_missing_project_for_all_projects_access_user(self):
        conn = IprConnection('user', 'pass')
        conn.get = mock.MagicMock(
            side_effect=[
                [{'name': 'OTHER'}],
                self._user_response(
                    groups=['create report access', 'all projects access'], projects=[]
                ),
            ]
        )
        conn.post = mock.MagicMock()

        conn.check_upload_permission('TEST')

        conn.post.assert_called_once_with('project', {'name': 'TEST'})
