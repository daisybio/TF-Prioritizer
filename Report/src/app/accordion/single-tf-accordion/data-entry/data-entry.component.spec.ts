import { ComponentFixture, TestBed } from '@angular/core/testing';

import { DataEntryComponent } from './data-entry.component';

describe('DataEntryComponent', () => {
  let component: DataEntryComponent;
  let fixture: ComponentFixture<DataEntryComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ DataEntryComponent ]
    })
    .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(DataEntryComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
