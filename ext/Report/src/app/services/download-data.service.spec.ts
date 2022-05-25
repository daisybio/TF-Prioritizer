import { TestBed } from '@angular/core/testing';

import { DownloadDataService } from './download-data.service';

describe('DownloadDataService', () => {
  let service: DownloadDataService;

  beforeEach(() => {
    TestBed.configureTestingModule({});
    service = TestBed.inject(DownloadDataService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });
});
