import { TestBed } from '@angular/core/testing';

import { OpenExternalUrlService } from './open-external-url.service';

describe('OpenExternalUrlService', () => {
  let service: OpenExternalUrlService;

  beforeEach(() => {
    TestBed.configureTestingModule({});
    service = TestBed.inject(OpenExternalUrlService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });
});
