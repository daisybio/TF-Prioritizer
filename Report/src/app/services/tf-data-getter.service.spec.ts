import { TestBed } from '@angular/core/testing';

import { TfDataGetterService } from './tf-data-getter.service';

describe('TfDataGetterService', () => {
  let service: TfDataGetterService;

  beforeEach(() => {
    TestBed.configureTestingModule({});
    service = TestBed.inject(TfDataGetterService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });
});
