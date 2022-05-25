import { TestBed } from '@angular/core/testing';

import { InformationGetterService } from './information-getter.service';

describe('InformationGetterService', () => {
  let service: InformationGetterService;

  beforeEach(() => {
    TestBed.configureTestingModule({});
    service = TestBed.inject(InformationGetterService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });
});
